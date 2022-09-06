#pragma once
#include <atomic>
#include <vector>
#include <Windows.Foundation.h>

class ParallelMesh;

// Our implementation for mesh decimation using a multi queue
class MultiQueue {
public:
	MultiQueue(int countVertices, int countQueues);
	void setErrors(std::vector<std::pair<float, int>> errorValues, int countThreads);

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
#endif // !SINGLE_THREADED_LOADING

	inline bool increaseAdjacentCollapses(int vertexId) {
		return 1 == ++m_adjacentToCollapses[vertexId];
	}

	inline void decreaseAdjacentCollapses(int vertexId) {
		m_adjacentToCollapses[vertexId]--;
		if(m_adjacentToCollapses[vertexId] < 0)
			__debugbreak();
	}

	inline bool getAdjacentCollapses(int vertexId) {
		//return true;
		return m_adjacentToCollapses[vertexId].load() == 1;
	}

	std::vector<int>& getNodes() {
		return m_nodes;
	}

	inline int getSize() {
		return static_cast<int>(m_nodes.size());
	}

private:
	int repairDown(int index, int queueId);
	bool repairUp(int index, int queueId);

	int getQueueId(int nodeId) {
		return nodeId / m_queueSize;
	}

	int getIndex(int nodeId, int queueId) {
		return nodeId - (queueId * m_queueSize);
	}

	int getNodeId(int index, int queueId) {
		return index + (queueId * m_queueSize);
	}
#ifndef SINGLE_THREADED_LOADING
	std::vector<CRITICAL_SECTION> m_vertexLocks; // for initialisation
#endif
	std::vector<CRITICAL_SECTION> m_queueLocks; // Lock per Queue

	int m_countQueues;
	int m_queueSize; // Amount of Vertices in one queue
	int m_size; // Amount of Vertices in all queues
	std::vector<std::atomic<float>> m_tmpError; // A temporary error for every vertex used to maintain the structure
	std::vector<int> m_nodes;
	std::vector<std::atomic<int>> m_nodeLookup; // A lookup to check which Vertex has which node
	ParallelMesh& m_mesh;
	std::vector<int> m_last; // Index of the Node that is to be deleted on pop, the last element in the queue;

	double m_summedError = 0;
	std::vector<std::atomic<int>> m_adjacentToCollapses; // increases by one if an collapse is performed on adjacent vertex or itself decreased by one when finished to prevent use in collapse when error is uncertain
};