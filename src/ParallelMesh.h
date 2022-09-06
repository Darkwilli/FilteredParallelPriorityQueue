#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <memory>
#include <filesystem>

#include "defines.h"
#include "Vertex.h"
#include "QuadricErrorMetric.h"

#include "PriorityQueueSingleThreaded.h"
#include "MultiQueue.h"
#include "FPPQ.h"
#include "FPPQMesh.h"

// Our mesh contains the vertices, the adjacency list for the face indices and surrounding vertices, the face index to corresponding vertex ids map and a pointer to a priority structure for mesh decimation
class ParallelMesh 
{
public:
	friend Vertex;
	friend PriorityQueueSingleThreaded;
	friend ErrorMetric;
	friend QuadricErrorMetric;
	friend MultiQueue;
	friend FPPQMesh;

	static ParallelMesh& getInstance() {
		static ParallelMesh instance;
		return instance;
	}
	// Loads the Mesh from the obj file in path
	bool loadFromObj(std::filesystem::path path);
	// Loads the mesh from exported data (because it is way faster then reloading the mesh from obj)
	void loadFromData(std::vector<Vertex>& vertices, std::vector<int>& adjacentVerticesIndices, std::vector<int>& faceIndices, std::vector<glm::ivec3>& faces);
	// Export the data so it can be loaded
	void exportData(std::vector<Vertex>& vertices, std::vector<int>& adjacentVerticesIndices, std::vector<int>& faceIndices, std::vector<glm::ivec3>& faces);

	// Mesh Decimation
	int reduceVerticesTo(float nrVertices, int countThreads);

	std::pair<std::vector<glm::vec3>, std::vector<glm::vec3>> getData();
	std::vector<glm::vec3> getVertices();
	std::vector<glm::ivec3> getFaces();

	// Initialise the QE Matrices
	void computeQuadricErrorMatrices(int countThreads);
	glm::vec3 getVertexPos(int index);

	Vertex& getVertex(int index);
	float getVertexError(int index);

	// Update the Priority Structure because the error from vertexId changed
	void updatePriorityStructure(int vertexId);

	void deletePriorityStructure() {
		m_priorityStructure = nullptr;
	}

	// Unlock the Vertex after the collapse
	void decreaseAdjacentCollapses(int vertexId)
	{
#ifndef SINGLE_THREADED
		m_priorityStructure->decreaseAdjacentCollapses(vertexId);
#else
		__debugbreak();
#endif // !SINGLE_THREADED
	}

	void debugCheckData(int activeException = -1);
	void reduceVertices(int finalSize, int countThreads); // Deletes the unused Vertices after decimation and updates the Vertex/Triangle Ids
private:
	void addTriangle(int currentIndex, int index1, int index2, int a, int b, int c, CRITICAL_SECTION& csF, CRITICAL_SECTION& csV, std::vector<std::pair<int, int>>& surplusVertices, std::vector<std::pair<int, glm::ivec3>>& surplusFaces);
	void addTriangleIndex(int triangleIndex, int currentIndex, int index1, int index2, int a, int b, int c, CRITICAL_SECTION& csF, CRITICAL_SECTION& csV, std::vector<std::pair<int, int>>& surplusVertices, std::vector<std::pair<int, int>>& surplusFaces);


	ParallelMesh() {}

	ParallelMesh(ParallelMesh const&) = delete;
	void operator=(ParallelMesh const&) = delete;

	// These are for the multi threaded compact operation (used for the reduceVertices operation after decimation)
	void initReplacementStart(int numThreads) 
	{
		m_replacementStart = std::vector<std::atomic<int>>(16 * numThreads); // Cache Optimisation Cache Line = 64B; because there is simultanious writing going on it costs us performance (2-5%) if the different entrys from m_replacementStart are in the same cache line.
		m_replacementEnd = std::vector<int>(numThreads);
	}
	void setReplacementStart(int newSize, int threadCount, int threadNr) 
	{
		int replaceMod = ( newSize) % threadCount;
		int replaceSlotSize, replaceStart;
		if (threadNr < replaceMod) {
			replaceSlotSize = (( newSize) / threadCount) + 1;
			replaceStart = 0 + replaceSlotSize * threadNr;
		}
		else {
			replaceSlotSize = (( newSize) / threadCount);
			replaceStart = 0 + replaceSlotSize * threadNr + replaceMod;
		}
		m_replacementStart[16 * threadNr] = replaceStart;
		m_replacementEnd[threadNr] = replaceStart + replaceSlotSize;
	}

	// The compact operation makes the m_vertices array dense again and creates a map for the changed vertex indices
	void compact(int oldSize, int newSize, int threadCount, int threadNr, std::vector<int>& vertexIndices, std::vector<int>& priorityNodes)
	{
		int mod = newSize % threadCount;
		int slotSize, start;
		if (threadNr < mod) {
			slotSize = (newSize / threadCount) + 1;
			start = slotSize * threadNr;
		}
		else {
			slotSize = (newSize / threadCount);
			start = slotSize * threadNr + mod;
		}
		int newVertex;
#ifndef SINGLE_THREADED
		int replacementPos;
#else
		int replacementPos = 0;
#endif
		int replacementIdx = threadNr;
		for (int i = start; i < start + slotSize; i++) {
			if (!m_vertices[i].isActive()) {
				newVertex = -1;
				while (newVertex == -1) {
#ifndef SINGLE_THREADED
					replacementPos = m_replacementStart[16 * replacementIdx]++;
					if (replacementPos >= m_replacementEnd[replacementIdx]) {
						replacementIdx = (replacementIdx + 1) % threadCount;
					}
					else 
#endif
					{
						if (priorityNodes[replacementPos] >= newSize)
							newVertex = priorityNodes[replacementPos];
					}
#ifdef SINGLE_THREADED
					replacementPos++;
#endif
				}
				m_vertices[i] = m_vertices[newVertex];
				vertexIndices[newVertex] = i;
			}
			else {
				vertexIndices[i] = i;
			}
			
		}
	}
	// Used for the compact operation
	std::vector<std::atomic<int>> m_replacementStart;
	std::vector<int> m_replacementEnd;

	// Mesh Information
	std::vector<Vertex> m_vertices;
	std::vector<int> m_adjacentVerticesIndices;
	std::vector<int> m_faceIndices;
	std::vector<glm::ivec3> m_faces;
	std::unique_ptr<PRIORITY_STRUCTURE> m_priorityStructure;
};

inline glm::vec3 ParallelMesh::getVertexPos(int index) {
	return m_vertices[index].getPos();
}

inline Vertex& ParallelMesh::getVertex(int index) {
	return m_vertices[index];
}

