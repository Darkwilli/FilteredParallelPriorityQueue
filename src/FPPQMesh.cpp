#include "FPPQMesh.h"
#include "Vertex.h"
#include "ParallelMesh.h"


FPPQMesh::FPPQMesh(int countVertices) : m_mesh(ParallelMesh::getInstance()), FPPQ<float>(countVertices)
{
	m_adjacentToCollapses = std::vector < std::atomic<int>>(m_size);
}

int FPPQMesh::pop()
{
	int popIndex = 0;
	int popIncreaseCount = maxPopIndex;
	int popLevel = 0;
	int popLevelBase = 0;
	int popLevelIndex = 0;
	int vertexId = -1;
	int lastIndex = m_last--; // We assume one does not pop the last vertex
	auto& adjacentVerticesIndices = m_mesh.m_adjacentVerticesIndices;
	while(vertexId == -1) {
		if(popIncreaseCount >= maxPopIndex || popIndex >= lastIndex) {
			//resets++;
			popIndex = 0;
			popLevel = 0;
			popIncreaseCount = 1;
			popLevelBase = 0;
			popLevelIndex = 0;
		} else {
			if(popLevelIndex >= (1 << (popLevel)) - 1) {
				popLevel++;
				popIncreaseCount++;
				popLevelIndex = 0;
				popLevelBase = (1 << (popLevel)) - 1;
				popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
			} else {
				popLevelIndex++;
				popIncreaseCount++;
				popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
			}

		}
		// Lock a root node
		while(vertexId == -1) {
			// If only Allowing root node replace tryLock loop with a lock
			while(!tryLock(popIndex))  // acquire lock
			{

				if(popIncreaseCount >= maxPopIndex || popIndex >= lastIndex) {
					popIndex = 0;
					popLevel = 0;
					popIncreaseCount = 1;
					popLevelBase = 0;
					popLevelIndex = 0;
				} else {
					if(popLevelIndex >= (1 << (popLevel)) - 1) {
						popLevel++;
						popIncreaseCount++;
						popLevelIndex = 0;
						popLevelBase = (1 << (popLevel)) - 1;
						popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
					} else {
						popLevelIndex++;
						popIncreaseCount++;
						popIndex = popLevelBase + ((popLevelIndex + omp_get_thread_num()) % (popLevelBase + 1));
					}

				}
			}
			vertexId = m_nodes[correctNodeId(popIndex)];
			if(m_nodeLookup[vertexId] == -1) {
				unlock(popIndex);
				__debugbreak();
				popIncreaseCount++;
				vertexId = -1;
			}
		}

		int lastValidPos = -1;

		if(increaseAdjacentCollapses(vertexId)) { // returns true if adjacent collapses == 0			
			{
				Vertex& vertex = m_mesh.getVertex(vertexId);
				// Prefetching increases the performance
				_mm_prefetch((char*)&vertex, _MM_HINT_T0);
				_mm_prefetch(((char*)&vertex) + 64, _MM_HINT_T0);
				_mm_prefetch(((char*)&vertex) + 128, _MM_HINT_T0);
				// returns true if adjacent collapses == 1 for surrounding vertices
				// (it is possible that 2 vertices with the same neighbour collapse at the same time resulting in an invalidation of the collapses) because of that the surrounding vertices must be checked
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
				if(possibleCollapse) {
					if(vertex.isValidCollapse())
					{
						m_nodeLookup[vertexId].store(-1);
						lock(lastIndex);
						int lastVertexId = m_nodes[correctNodeId(lastIndex)];
						m_nodeLookup[lastVertexId].store(popIndex);
						m_nodes[correctNodeId(popIndex)] = (lastVertexId);
						m_notRepairedUp[lastVertexId] += 1;
						m_nodes[correctNodeId(lastIndex)] = (-1);
						unlock(lastIndex);
						int nextNodeId = popIndex;
						while(nextNodeId != -1) {
							nextNodeId = repairDown(nextNodeId);
						}
						m_notRepairedUp[lastVertexId] -= 1;
					}
					else {
#ifdef NEXT_BEST_ERROR
						m_notRepairedUp[vertexId] += 1;
						float newError = vertex.setNextBestError();
						m_tmpError[vertexId] = newError;

#else
						if(m_tmpError[vertexId] >= 1000000000000000000000000.0)
							__debugbreak();
						m_notRepairedUp[vertexId] += 1;
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
						int nextNodeId = popIndex;
						while(nextNodeId != -1) {
							nextNodeId = repairDown(nextNodeId);
						}
						m_notRepairedUp[vertexId] -= 1;
						vertexId = -1;
						popIncreaseCount--;
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
					unlock(popIndex);
				}
			}
		} else {
			decreaseAdjacentCollapses(vertexId);
			vertexId = -1;
			unlock(popIndex);
		}
		popIncreaseCount++;
	}
	return vertexId;

}

void FPPQMesh::update(int vertexId, bool fromInsertion) {
	float errorValue = 0.f;
	if(!fromInsertion)
		errorValue = m_mesh.getVertexError(vertexId);
	FPPQ<float>::update(vertexId, errorValue, fromInsertion);
}
