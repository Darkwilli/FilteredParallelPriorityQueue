#pragma once

#include "FPPQ.h"

class ParallelMesh;

class FPPQMesh : public FPPQ<float> {
public:
	FPPQMesh(int countVertices);
	virtual int pop() override;

	void setErrors(std::vector<std::pair<float, int>> errorValues, int countThreads) {
		FPPQ<float>::fill(errorValues, countThreads);
	}

	virtual void update(int vertexId, bool fromInsertion = false);

	inline bool increaseAdjacentCollapses(int vertexId) {
		return 1 == ++m_adjacentToCollapses[vertexId];
	}

	inline void decreaseAdjacentCollapses(int vertexId) {
		m_adjacentToCollapses[vertexId]--;
		if(m_adjacentToCollapses[vertexId] < 0)
			__debugbreak();
	}

	inline bool getAdjacentCollapses(int vertexId) {
		return m_adjacentToCollapses[vertexId].load() == 1;
	}

private:
	// Stored them here because they are needed here most of the time
	std::vector<std::atomic<int>> m_adjacentToCollapses; // increases by one if a collapse is performed on an adjacent vertex or itself. Decreased by one when finished. To prevent use of vertex while collapsing
	ParallelMesh& m_mesh;

};