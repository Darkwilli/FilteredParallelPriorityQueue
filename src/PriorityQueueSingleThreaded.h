#pragma once
#include <vector>
#include <Windows.Foundation.h>
class ParallelMesh;

// A fast single threaded priority queue implementation for mesh decimation
class PriorityQueueSingleThreaded {
public:
	PriorityQueueSingleThreaded(int countVertices);
	void setErrors(std::vector<std::pair<float, int>> errorValues);

	int pop(); // returns vertexId of node to pop
	void update(int vertexId);

	void debugCheckHeap(); // Checks if the Heap is still valid
#ifndef SINGLE_THREADED_LOADING
	void lock(int vertexId) {
		EnterCriticalSection(&m_vertexLocks[vertexId]);
	}

	void unlock(int vertexId) {
		LeaveCriticalSection(&m_vertexLocks[vertexId]);
	}
#endif
	std::vector<int>& getNodes() {
		return m_nodes;
	}

private:
	int repairDown(int nodeId);
	bool repairUp(int nodeId);
#ifndef SINGLE_THREADED_LOADING
	std::vector<CRITICAL_SECTION> m_vertexLocks; // for initialisation
#endif
	int m_size;
	std::vector<float> m_tmpError; // A temporary error for every vertex used to maintain the structure
	std::vector<int> m_nodes; // Left child: m_nodeId * 2 + 1; Right child: m_nodeId * 2 + 2; Parent: (m_nodeId - 1) / 2 ;
	std::vector<int> m_nodeLookup; // A lookup to check which Vertex has which node
	ParallelMesh& m_mesh;
	int m_last; // Index of the Node that is to be deleted on pop, the last element in the array;
	int m_repairCount;
	int m_repairDirection;

	double m_summedError = 0;
};