#include <algorithm>
#include <cassert>

#include "PriorityQueueSingleThreaded.h"
#include "Vertex.h"
#include "ParallelMesh.h"

#define PARENT_NODE_ID ((nodeId - 1 ) / 2)
#define LEFT_CHILD_NODE_ID (nodeId * 2 + 1)
#define RIGHT_CHILD_NODE_ID (nodeId * 2 + 2)

PriorityQueueSingleThreaded::PriorityQueueSingleThreaded(int countVertices) : m_mesh(ParallelMesh::getInstance())
{
	m_size = countVertices;
	m_nodes = std::vector<int>(m_size);
	m_tmpError = std::vector<float>(m_size);

	m_nodeLookup = std::vector<int>(m_size);

	m_last = m_size - 1;
	m_repairCount = 0;
#ifndef SINGLE_THREADED_LOADING
	m_vertexLocks = std::vector<CRITICAL_SECTION>(m_size);
	int spincount = 0x01001000;
	for (int i = 0; i < m_size; i++) {
		if (!InitializeCriticalSectionAndSpinCount(&m_vertexLocks[i], spincount))
			__debugbreak();
	}
#endif

}

void PriorityQueueSingleThreaded::setErrors(std::vector<std::pair<float, int>> errorValues)
{

	std::sort(errorValues.begin(), errorValues.end(), std::less<std::pair<float, int>>());
	for (int i = 0; i < m_size; i++) {
		m_nodes[i] = errorValues[i].second;
		m_tmpError[errorValues[i].second] = errorValues[i].first;
		m_nodeLookup[errorValues[i].second] = i;
	}
}
int PriorityQueueSingleThreaded::pop()
{
	int vertexId = -1;
	while (vertexId == -1)
	{
		vertexId = m_nodes[0];
		if (!m_mesh.getVertex(vertexId).isValidCollapse())
		{
#ifdef NEXT_BEST_ERROR
			float newError = m_mesh.getVertex(vertexId).setNextBestError();
			m_tmpError[vertexId] = newError;
#else
			if (m_tmpError[vertexId] >= 1000000000000000000000000.0)
				__debugbreak();
			m_tmpError[vertexId] +=  10000000000000000000000.0;
#endif
			int nextNodeId = 0;
			while (nextNodeId != -1) {
				nextNodeId = repairDown(nextNodeId);
			}
			vertexId = -1;
		}
		
	}	

	m_nodeLookup[vertexId] = -1;

	int lastVertexId = m_nodes[m_last];
	m_nodes[m_last] = -1;
	m_last--;
	m_nodeLookup[lastVertexId] = 0;
	m_nodes[0] = lastVertexId;

	int nextNodeId = 0;
	while (nextNodeId != -1) {
		nextNodeId = repairDown(nextNodeId);
	}
	return vertexId;
}

void PriorityQueueSingleThreaded::update(volatile int vertexId)
{
	float trueError = m_mesh.getVertexError(vertexId);
	//m_mesh->getVertex(vertexId).unlock();
	bool directionUp = false;
	if (trueError < m_tmpError[vertexId]) { // We now have a smaller Error (We need to look upwards)
		directionUp = true;
	}
	m_tmpError[vertexId] = trueError;

	volatile int nodeId = m_nodeLookup[vertexId];
	if (nodeId < 0)
		__debugbreak();
	if (nodeId > m_last)
		__debugbreak();
	if (directionUp) {
		bool finished = false;
		while (!finished)
		{
			finished = repairUp(nodeId);
			nodeId = PARENT_NODE_ID;
			/*if ((m_nodeLookup[vertexId] != nodeId) && !finished)
				__debugbreak();*/
		}
		//m_locks[nodeId].clear(std::memory_order_release);
	}
	else { // downwards
		int nextNodeId = m_nodeLookup[vertexId];
		while (nextNodeId != -1) {
			nextNodeId = repairDown(nextNodeId);
		}
	}
}

void PriorityQueueSingleThreaded::debugCheckHeap()
{
	volatile int nodeId, parentNodeId;
	for (nodeId = 1; nodeId < m_nodes.size(); nodeId++) {
		if (nodeId > m_last - 2) // -2 because we ignore the next childs when they can be removed in the next 2 pops
			return;

		parentNodeId = PARENT_NODE_ID;

		if (m_tmpError[m_nodes[nodeId]] < m_tmpError[m_nodes[parentNodeId]])
		{
			volatile int vidc = m_nodes[nodeId];
			volatile int vidp = m_nodes[parentNodeId];
			volatile float errtc = m_tmpError[vidc];
			volatile float errtp = m_tmpError[vidp];
			//float errdbc = m_nodes[i].getDebugError();
			//float errdbp = m_nodes[(i - 1) / 2].getDebugError();
			volatile float errvc = m_mesh.getVertexError(m_nodes[nodeId]);
			volatile float errvp = m_mesh.getVertexError(m_nodes[parentNodeId]);
			__debugbreak();
		}
		if (m_tmpError[m_nodes[nodeId]] != m_mesh.getVertexError(m_nodes[nodeId]) && m_tmpError[m_nodes[nodeId]] < 100000000000000000000.0)
		{
			volatile int vidc = m_nodes[nodeId];
			volatile int vidp = m_nodes[parentNodeId];
			volatile float errtc = m_tmpError[vidc];
			volatile float errtp = m_tmpError[vidp];
			//float errdbc = m_nodes[i].getDebugError();
			//float errdbp = m_nodes[(i - 1) / 2].getDebugError();
			volatile float errvc = m_mesh.getVertexError(m_nodes[nodeId]);
			volatile float errvp = m_mesh.getVertexError(m_nodes[parentNodeId]);
			__debugbreak();
		}
	}
}

int PriorityQueueSingleThreaded::repairDown(int nodeId)
{
	volatile int leftChildNodeId, rightChildNodeId;

	leftChildNodeId = LEFT_CHILD_NODE_ID;
	rightChildNodeId = RIGHT_CHILD_NODE_ID;	
	volatile int nextNodeId;
	if (rightChildNodeId <= m_last) {
		volatile int vertexIdLeft = m_nodes[leftChildNodeId];
		volatile int vertexIdRight = m_nodes[rightChildNodeId];

		volatile float valueLeft = m_tmpError[vertexIdLeft];
		volatile float valueRight = m_tmpError[vertexIdRight];
		volatile int vertexId = m_nodes[nodeId];
		if (valueLeft >= valueRight) {
			if (valueRight < m_tmpError[vertexId]){
				nextNodeId = rightChildNodeId;
				int vertexIdChild = vertexIdRight;
				m_nodeLookup[vertexId] = nextNodeId;
				m_nodes[nextNodeId] = (vertexId);
				m_nodeLookup[vertexIdChild] = nodeId;
				m_nodes[nodeId] = vertexIdChild;
			}
			else {
				return -1;
			}
		}
		else if (valueLeft < m_tmpError[vertexId]) {
			nextNodeId = leftChildNodeId;
			volatile int vertexIdChild = vertexIdLeft;
			m_nodeLookup[vertexId] = nextNodeId;
			m_nodes[nextNodeId] = (vertexId);
			m_nodeLookup[vertexIdChild] = nodeId;
			m_nodes[nodeId] = vertexIdChild;
		} else {
			return -1;
		}
	}
	else {
		if (leftChildNodeId == m_last) {
			volatile int vertexIdLeft = m_nodes[leftChildNodeId];
			volatile float valueLeft = m_tmpError[vertexIdLeft];
			volatile int vertexId = m_nodes[nodeId];
			if (valueLeft < m_tmpError[vertexId]) {
				nextNodeId = leftChildNodeId;
				volatile int vertexIdChild = vertexIdLeft;
				m_nodeLookup[vertexId] = nextNodeId;
				m_nodes[nextNodeId] = (vertexId);
				m_nodeLookup[vertexIdChild] = nodeId;
				m_nodes[nodeId] = vertexIdChild;
			}
			else {
				return -1;
			}
		}
		else {
			return -1;
		}
	}
	return nextNodeId;
}

bool PriorityQueueSingleThreaded::repairUp(int nodeId)
{
	int parentNodeId = PARENT_NODE_ID;
	int parentVertexId = m_nodes[parentNodeId];
	float valueParent = m_tmpError[parentVertexId];
	int vertexId = m_nodes[nodeId];
	if (m_tmpError[vertexId] < valueParent) { // We are now smaller than the parent

		m_nodeLookup[vertexId] = parentNodeId; // swap
		m_nodes[parentNodeId] = vertexId;
		m_nodeLookup[parentVertexId] = nodeId;
		m_nodes[nodeId] = parentVertexId;

		if (parentNodeId == 0)
			return true; // The node is now the root so we are finished
		return false; // We have to look at our parent
	}
	else {
		// we have the correct constellation
		return true;
	}
}
