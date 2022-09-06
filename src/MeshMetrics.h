#pragma once

#include <vector>
#include <unordered_map>
#include <omp.h>
#include <atomic>
#include <glm/ext.hpp>



namespace MeshMetrics {
	void computeHausdorff(const std::vector<glm::vec3> verticesA, const std::vector<glm::ivec3> facesA, const std::vector<glm::vec3> verticesB, const std::vector<glm::ivec3> facesB, double& mean, double& max)
	{
		// We used a library for Hausdorff calculation which is not included
		
		mean = -1;
		max = -1;
		return;
	}


	struct ivec3Hash {
	public:
		size_t operator()(const glm::ivec3 p) const
		{
			return std::hash<int>()(p.x) ^ std::hash<int>()(p.y) ^ std::hash<int>()(p.z);
		}
	};

	struct vec3Hash {
	public:
		size_t operator()(const glm::vec3 p) const
		{
			return std::hash<float>()(p.x) ^ std::hash<float>()(p.y) ^ std::hash<float>()(p.z);
		}
	};

	// Compute the ratio of same Vertices/ Faces for two decimated Meshes
	void computeSimilarity(const std::vector<glm::vec3> verticesA, const std::vector<glm::ivec3> facesA, const std::vector<glm::vec3> verticesB, const std::vector<glm::ivec3> facesB, double& vRatio, double& fRatio) {

		std::vector<int> mapBtoA(verticesB.size(), -1);

		std::unordered_map<glm::vec3, int, vec3Hash> aVerts;

		std::vector<glm::vec3> bCopy = verticesB;

		for (int a = 0; a < verticesA.size(); a++) {
			aVerts.insert({ verticesA[a], a });
		}

		std::atomic<int> sameVertices = 0;
//#pragma omp parallel for num_threads(32)
		for (int b = 0; b < verticesB.size(); b++) {
			if (aVerts.contains(verticesB[b])) {
				mapBtoA[b] = aVerts[verticesB[b]];
				sameVertices++;
				//break;
			}
		}

		vRatio = double(sameVertices) / verticesA.size();

		std::unordered_map<glm::ivec3, int, ivec3Hash> aFaces;

		for (int a = 0; a < facesA.size(); a++) {
			aFaces.insert({ facesA[a], a });
		}

		std::atomic<int> sameFaces = 0;
//#pragma omp parallel for num_threads(32) schedule(guided)		
		for(int b = 0; b < facesB.size(); b++) {
			if(aFaces.contains(glm::ivec3(mapBtoA[facesB[b].x], mapBtoA[facesB[b].y], mapBtoA[facesB[b].z]))) 
			{
				sameFaces++;
				//break;
			}
		}
		

		fRatio = double(sameFaces) / facesA.size();
	}

}
