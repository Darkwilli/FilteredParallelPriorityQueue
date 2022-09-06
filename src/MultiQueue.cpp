#include <algorithm>
#include <execution>
#include <cassert>
#include <random>

#include "MultiQueue.h"
#include "Vertex.h"
#include "ParallelMesh.h"


#define PARENT_NODE_ID ((nodeId - 1 ) / 2)
#define LEFT_CHILD_NODE_ID (nodeId * 2 + 1)
#define RIGHT_CHILD_NODE_ID (nodeId * 2 + 2)

MultiQueue::MultiQueue(int countVertices, int countQueues) : m_mesh(ParallelMesh::getInstance()) {
	m_countQueues = countQueues;
	m_queueSize = countVertices / m_countQueues;
	int queueMod = countVertices % m_countQueues;
	queueMod == 0 ? m_queueSize : m_queueSize++;

	m_size = m_queueSize * m_countQueues;
	m_nodes = std::vector<int>(m_size, -1);
	m_tmpError = std::vector<std::atomic<float>>(m_size);

	m_nodeLookup = std::vector<std::atomic<int>>(m_size);
	m_adjacentToCollapses = std::vector < std::atomic<int>>(m_size);
	m_last = std::vector<int>(m_countQueues);

	for(int i = 0; i < m_countQueues; i++) {
		int additionalElement;
		if(queueMod != 0)
			additionalElement = queueMod > i ? 0 : 1;
		else
			additionalElement = 0;
		m_last[i] = m_queueSize * (i + 1) - additionalElement - 1;
	}
	int spincount = 0x01001000;
#ifndef SINGLE_THREADED_LOADING
	m_vertexLocks = std::vector<CRITICAL_SECTION>(m_size);
	for(int i = 0; i < m_size; i++) {
		if(!InitializeCriticalSectionAndSpinCount(&m_vertexLocks[i], spincount))
			__debugbreak();
	}
#endif
	m_queueLocks = std::vector<CRITICAL_SECTION>(m_countQueues);
	for(int i = 0; i < m_countQueues; i++) {
		if(!InitializeCriticalSectionAndSpinCount(&m_queueLocks[i], spincount))
			__debugbreak();
	}

}

void MultiQueue::setErrors(std::vector<std::pair<float, int>> errorValues, int countThreads) {
	std::sort(std::execution::par_unseq, errorValues.begin(), errorValues.end(), std::less<std::pair<float, int>>());
#pragma omp parallel for num_threads(countThreads)
	for(int i = 0; i < errorValues.size(); i++) {
		int queueIndex = i / m_countQueues;
		int queueId = i % m_countQueues;

		int index = queueId * m_queueSize + queueIndex;


		m_nodes[index] = errorValues[i].second;
		m_tmpError[errorValues[i].second] = errorValues[i].first;
		m_nodeLookup[errorValues[i].second] = index;
	}
}
int MultiQueue::pop() {

	thread_local static std::random_device rd;
	//thread_local static PCG64 rng(rd());
	thread_local static std::mt19937 rng(rd());
	//thread_local static std::mt19937 rng(0);
	//thread_local static std::mt19937 rng(omp_get_thread_num());
	std::uniform_int_distribution<int> dist;
	dist = std::uniform_int_distribution<int>(0, m_countQueues - 1);

	int popIndex = 0;
	int vertexId = -1;
	auto& adjacentVerticesIndices = m_mesh.m_adjacentVerticesIndices;
	while(vertexId == -1) {
		// Lock two random queues
		int queueId1 = dist(rng);
		while(!TryEnterCriticalSection(&m_queueLocks[queueId1]))  // acquire lock
		{
			queueId1 = dist(rng);
		}
		int queueId2 = dist(rng);
		while(queueId2 == queueId1 || !TryEnterCriticalSection(&m_queueLocks[queueId2]))  // acquire lock
		{
			queueId2 = dist(rng);
		}
		
		// Get the better of the two root elements

		int vertexId1 = m_nodes[queueId1 * m_queueSize];
		int vertexId2 = m_nodes[queueId2 * m_queueSize];

		float err1 = m_tmpError[vertexId1];
		float err2 = m_tmpError[vertexId2];

		int queueId;
		if(err1 < err2) {
			vertexId = vertexId1;
			LeaveCriticalSection(&m_queueLocks[queueId2]);
			queueId = queueId1;
		} else {
			vertexId = vertexId2;
			LeaveCriticalSection(&m_queueLocks[queueId1]);
			queueId = queueId2;
		}
		popIndex = queueId * m_queueSize;
		int lastValidPos = -1;
		if(increaseAdjacentCollapses(vertexId)) { // returns true if adjacent collapses == 0			
			{
				Vertex& vertex = m_mesh.getVertex(vertexId);
				_mm_prefetch((char*)&vertex, _MM_HINT_T0);
				_mm_prefetch(((char*)&vertex) + 64, _MM_HINT_T0);
				_mm_prefetch(((char*)&vertex) + 128, _MM_HINT_T0);
				bool possibleCollapse = true;
				{
					int startPos = (slotSize + 1) * vertexId;
					int start = adjacentVerticesIndices[startPos];
					int currentPos = startPos + 1;
					int ctr = 0;
					bool stop = false;
					while(!stop && adjacentVerticesIndices[currentPos] != -1) {
						if(adjacentVerticesIndices[currentPos] != vertexId)
							if(!increaseAdjacentCollapses(adjacentVerticesIndices[currentPos])) // increase to mark them used (decrease will be after the collapse or in case someone else was faster)
							{
								lastValidPos = currentPos;
								possibleCollapse = false;
								break;
							}
						if(++ctr != slotSize) {
							currentPos += 1;
						} else {
							ctr = 0;
							if(start != -1) {
								startPos = start;
								start = adjacentVerticesIndices[startPos];
								currentPos = startPos + 1;
							} else {
								stop = true;
							}
						}
					}
				}
				// returns true if adjacent collapses == 1 for itself and vertex to collapse to 
				// (it is possible that 2 vertices with the same neighbour collapse at the same time resulting in an invalidation of the collapses) because of that the vertex to collapse to must be chacked as well
				if(possibleCollapse) {
					if(vertex.isValidCollapse())
					{
						m_nodeLookup[vertexId].store(-1);
						int lastVertexId = m_nodes[m_last[queueId]];
						m_nodeLookup[lastVertexId].store(popIndex);
						m_nodes[popIndex] = (lastVertexId);
						m_nodes[m_last[queueId]] = (-1);
						m_last[queueId]--;

						int nextNodeId = 0;
						while(nextNodeId != -1) {
							nextNodeId = repairDown(nextNodeId, queueId);
						}

					}
					else {
#ifdef NEXT_BEST_ERROR
						float newError = vertex.setNextBestError();
						m_tmpError[vertexId] = newError;

#else
						if(m_tmpError[vertexId] >= 1000000000000000000000000.0)
							__debugbreak();
						m_tmpError[vertexId] += 10000000000000000000000.0;

#endif // NEXT_BEST_ERROR

						int startPos = (slotSize + 1) * vertexId;
						int start = adjacentVerticesIndices[startPos];
						int currentPos = startPos + 1;
						int ctr = 0;
						bool stop = false;
						while(!stop && adjacentVerticesIndices[currentPos] != -1) {
							decreaseAdjacentCollapses(adjacentVerticesIndices[currentPos]); // Decrease this vertex will not be collapsed
							if(++ctr != slotSize) {
								currentPos += 1;
							} else {
								ctr = 0;
								if(start != -1) {
									startPos = start;
									start = adjacentVerticesIndices[startPos];
									currentPos = startPos + 1;
								} else {
									stop = true;
								}
							}
						}
						int nextNodeId = 0;
						while(nextNodeId != -1) {
							nextNodeId = repairDown(nextNodeId, queueId);
						}
						vertexId = -1;
					}
				} else {
					int startPos = (slotSize + 1) * vertexId;
					int start = adjacentVerticesIndices[startPos];
					int currentPos = startPos + 1;
					int ctr = 0;
					bool stop = false;
					while(!stop && adjacentVerticesIndices[currentPos] != -1) {
						decreaseAdjacentCollapses(adjacentVerticesIndices[currentPos]); // Decrease this vertex will not be collapsed
						if(currentPos == lastValidPos) {
							break;
						}
						if(++ctr != slotSize) {
							currentPos += 1;
						} else {
							ctr = 0;
							if(start != -1) {
								startPos = start;
								start = adjacentVerticesIndices[startPos];
								currentPos = startPos + 1;
							} else {
								stop = true;
							}
						}
					}
					vertexId = -1;
				}
			}
		} else {
			decreaseAdjacentCollapses(vertexId);
			vertexId = -1;
		}
		LeaveCriticalSection(&m_queueLocks[queueId]);
	}

	return vertexId;
}

void MultiQueue::update(volatile int vertexId) {
	
	int vertexNodeId = m_nodeLookup[vertexId];
	int queueId = getQueueId(vertexNodeId);

	EnterCriticalSection(&m_queueLocks[queueId]);
	vertexNodeId = m_nodeLookup[vertexId];
	float trueError = m_mesh.getVertexError(vertexId);
	//m_mesh->getVertex(vertexId).unlock();
	bool directionUp = false;
	if(trueError < m_tmpError[vertexId]) { // We now have a smaller Error (We need to look upwards)
		directionUp = true;
	}
	m_tmpError[vertexId] = trueError;

	int nodeId = getIndex(vertexNodeId, queueId);
	if(nodeId < 0)
		__debugbreak();
	if(directionUp) {
		bool finished = false;
		while(!finished) {
			finished = repairUp(nodeId, queueId);
			nodeId = PARENT_NODE_ID;
			/*if ((m_nodeLookup[vertexId] != nodeId) && !finished)
				__debugbreak();*/
		}
		//m_locks[nodeId].clear(std::memory_order_release);
	} else { // downwards
		int nextNodeId = nodeId;
		while(nextNodeId != -1) {
			nextNodeId = repairDown(nextNodeId, queueId);
		}
	}
	LeaveCriticalSection(&m_queueLocks[queueId]);
	//debugCheckHeap();
}

void MultiQueue::debugCheckHeap() {
	// We dont do that here
}

int MultiQueue::repairDown(int nodeId, int queueId) {
	volatile int leftChildNodeId, rightChildNodeId;

	int offset = queueId * m_queueSize;

	leftChildNodeId = LEFT_CHILD_NODE_ID;
	rightChildNodeId = RIGHT_CHILD_NODE_ID;
	volatile int nextNodeId;
	if(rightChildNodeId + offset <= m_last[queueId]) {
		volatile int vertexIdLeft = m_nodes[leftChildNodeId + offset];
		volatile int vertexIdRight = m_nodes[rightChildNodeId + offset];

		volatile float valueLeft = m_tmpError[vertexIdLeft];
		volatile float valueRight = m_tmpError[vertexIdRight];
		volatile int vertexId = m_nodes[nodeId + offset];
		if(valueLeft >= valueRight) {
			if(valueRight < m_tmpError[vertexId]) {
				nextNodeId = rightChildNodeId;
				int vertexIdChild = vertexIdRight;
				m_nodeLookup[vertexId] = nextNodeId + offset;
				m_nodes[nextNodeId + offset] = (vertexId);
				m_nodeLookup[vertexIdChild] = nodeId + offset;
				m_nodes[nodeId + offset] = vertexIdChild;
			} else {
				return -1;
			}
		} else if(valueLeft < m_tmpError[vertexId]) {
			nextNodeId = leftChildNodeId;
			volatile int vertexIdChild = vertexIdLeft;
			m_nodeLookup[vertexId] = nextNodeId + offset;
			m_nodes[nextNodeId + offset] = (vertexId);
			m_nodeLookup[vertexIdChild] = nodeId + offset;
			m_nodes[nodeId + offset] = vertexIdChild;
		} else {
			return -1;
		}
	} else {
		if(leftChildNodeId + offset == m_last[queueId]) {
			volatile int vertexIdLeft = m_nodes[leftChildNodeId + offset];
			volatile float valueLeft = m_tmpError[vertexIdLeft];
			volatile int vertexId = m_nodes[nodeId + offset];
			if(valueLeft < m_tmpError[vertexId]) {
				nextNodeId = leftChildNodeId;
				volatile int vertexIdChild = vertexIdLeft;
				m_nodeLookup[vertexId] = nextNodeId + offset;
				m_nodes[nextNodeId + offset] = (vertexId);
				m_nodeLookup[vertexIdChild] = nodeId + offset;
				m_nodes[nodeId + offset] = vertexIdChild;
			} else {
				return -1;
			}
		} else {
			return -1;
		}
	}
	return nextNodeId;
}

bool MultiQueue::repairUp(int nodeId, int queueId) {
	int offset = queueId * m_queueSize;
	int parentNodeId = PARENT_NODE_ID;
	int parentVertexId = m_nodes[parentNodeId + offset];
	float valueParent = m_tmpError[parentVertexId];
	int vertexId = m_nodes[nodeId + offset];
	if(m_tmpError[vertexId] < valueParent) { // We are now smaller than the parent

		m_nodeLookup[vertexId] = parentNodeId + offset; // swap
		m_nodes[parentNodeId + offset] = vertexId;
		m_nodeLookup[parentVertexId] = nodeId + offset;
		m_nodes[nodeId + offset] = parentVertexId;

		if(parentNodeId == 0)
			return true; // The node is now the root so we are finished
		return false; // We have to look at our parent
	} else {
		// we have the correct constellation
		return true;
	}
}
