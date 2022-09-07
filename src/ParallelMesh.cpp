
#include <unordered_map>
#include <cassert>
#include <atomic>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>

#include <execution>
#include <algorithm>

#include "ParallelMesh.h"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader.h"

#undef min

// Reduce the mesh to nrVertices (if it is smaller than 1 it is treated as a percentage) and returns the amount of leftover Vertices
int ParallelMesh::reduceVerticesTo(float nrVertices, int countThreads) {
	return reduceVerticesTo(static_cast<int>(m_vertices.size() * nrVertices), countThreads);
}
int ParallelMesh::reduceVerticesTo(int nrVertices, int countThreads)
{
	if (nrVertices >= m_vertices.size())
		return nrVertices;
	int finalVertexCount = nrVertices;

#ifndef SINGLE_THREADED
	int num_threads = countThreads;
	std::atomic<int> num_cur_threads = num_threads;
#else
	int num_threads = 1;
	int num_cur_threads = num_threads;
	int removedVertices = 0;
#endif
	int size = static_cast<int>(m_vertices.size());

#pragma omp parallel num_threads(num_threads)
	{
#pragma omp for nowait schedule(guided)
		for (int i = 0; i < m_vertices.size() - finalVertexCount; i++)
		{	
			int collapseVertexIndex = m_priorityStructure->pop();
			m_vertices[collapseVertexIndex].collapse(collapseVertexIndex);
		}
		num_cur_threads--;
#pragma omp barrier
	}
#ifndef BENCHMARK
	m_priorityStructure->debugCheckHeap();
#endif // !BENCHMARK

	return finalVertexCount;
}

bool ParallelMesh::loadFromObj(std::filesystem::path path) {
	if(!std::filesystem::exists(path)) {
		return false;
	}
	std::string inputfile = path.string();
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());

	if (!err.empty()) { // `err` may contain warning message.
		std::cerr << err << std::endl;
	}

	if (!ret) {
		return false;
	}

	auto& vertices = attrib.GetVertices();

	m_vertices.resize(vertices.size() / 3);

#ifdef MULTI_QUEUE
	m_priorityStructure = std::make_unique<PRIORITY_STRUCTURE>(static_cast<int>(vertices.size()) / 3, NUM_THREADS * 2);
#else
	m_priorityStructure = std::make_unique<PRIORITY_STRUCTURE>(static_cast<int>(vertices.size()) / 3);
#endif // MULTI_QUEUE



	CRITICAL_SECTION csF;
	if (!InitializeCriticalSectionAndSpinCount(&csF,
		0x01001000))
		return false;
	CRITICAL_SECTION csV;
	if (!InitializeCriticalSectionAndSpinCount(&csV,
		0x01001000))
		return false;
	std::vector<std::pair<int, int>> surplusVertices;
	int n = (static_cast<int>(vertices.size()) / 3) * (slotSize + 1);
	m_adjacentVerticesIndices = std::vector<int>(n, -1); // 1 continue 10 slots
	std::vector<std::pair<int, int>> surplusFaces;
	m_faceIndices = std::vector<int>(n, -1); // 1 continue 10 slots
	m_faces = std::vector<glm::ivec3>(shapes[0].mesh.indices.size() / 3); // 1 continue 10 slots

#ifndef SINGLE_THREADED_LOADING
#pragma omp parallel num_threads(32)
#endif // !SINGLE_THREADED_LOADING
	{
#ifndef SINGLE_THREADED_LOADING
#pragma omp for schedule(guided)
#endif // !SINGLE_THREADED_LOADING
		for (int v = 0; v < vertices.size() / 3; v++) {
			m_vertices[v] = (Vertex({ (float)vertices[3 * v + 0] , (float)vertices[3 * v + 1], (float)vertices[3 * v + 2] }, v));
		}

		if (shapes.size() != 1)
			return false;
#ifndef SINGLE_THREADED_LOADING
#pragma omp barrier
#pragma omp for schedule(guided)
#endif // !SINGLE_THREADED_LOADING
		for (int f = 0; f < shapes[0].mesh.indices.size() / 3; f++) {
			int a = (int)shapes[0].mesh.indices[3 * f].vertex_index;
			int b = (int)shapes[0].mesh.indices[3 * f + 1].vertex_index;
			int c = (int)shapes[0].mesh.indices[3 * f + 2].vertex_index;

#ifndef SINGLE_THREADED_LOADING
			m_priorityStructure->lock(a);
			addTriangleIndex(f, a, b, c, a, b, c, csF, csV, surplusVertices, surplusFaces);
			m_priorityStructure->unlock(a);
			m_priorityStructure->lock(b);
			addTriangleIndex(f, b, a, c, a, b, c, csF, csV, surplusVertices, surplusFaces);
			m_priorityStructure->unlock(b);
			m_priorityStructure->lock(c);
			addTriangleIndex(f, c, a, b, a, b, c, csF, csV, surplusVertices, surplusFaces);
			m_priorityStructure->unlock(c);
#else
			addTriangleIndex(f, a, b, c, a, b, c, csF, csV, surplusVertices, surplusFaces);
			addTriangleIndex(f, b, a, c, a, b, c, csF, csV, surplusVertices, surplusFaces);
			addTriangleIndex(f, c, a, b, a, b, c, csF, csV, surplusVertices, surplusFaces);
#endif
			m_faces[f] = { a,b,c };
		}
#ifndef SINGLE_THREADED_LOADING
#pragma omp single nowait
#endif // !SINGLE_THREADED_LOADING
		{
			for(auto& vertex : surplusVertices) {
				bool finished = false;
				int startPos = (slotSize + 1) * vertex.first;
				int start = m_adjacentVerticesIndices[startPos];
				bool hasAdditionalSpace = false;
				if(start != -1) {
					startPos = start;
					start = m_adjacentVerticesIndices[start];
					int currentPos = startPos + 1;
					bool alreadyInserted = false;
					int ctr = 0;
					while(!alreadyInserted && m_adjacentVerticesIndices[currentPos] != -1) {
						if(m_adjacentVerticesIndices[currentPos] == vertex.second) {
							alreadyInserted = true;
						}
						if(++ctr != slotSize) {
							currentPos += 1;
						} else {
							ctr = 0;
							if(start != -1) {
								startPos = start;
								start = m_adjacentVerticesIndices[startPos];
								currentPos = startPos + 1;
							} else {
								if(!alreadyInserted) {
									m_adjacentVerticesIndices[startPos] = static_cast<int>(m_adjacentVerticesIndices.size());
									m_adjacentVerticesIndices.push_back(-1);
									m_adjacentVerticesIndices.push_back(vertex.second);
									for(int i = 0; i < slotSize - 1; i++) {
										m_adjacentVerticesIndices.push_back(-1);
									}
									alreadyInserted = true;
								}
							}
						}
					}
					if(!alreadyInserted) {
						m_adjacentVerticesIndices[currentPos] = vertex.second;
					}
				} else {
					m_adjacentVerticesIndices[startPos] = static_cast<int>(m_adjacentVerticesIndices.size());
					m_adjacentVerticesIndices.push_back(-1);
					m_adjacentVerticesIndices.push_back(vertex.second);
					for(int i = 0; i < slotSize - 1; i++) {
						m_adjacentVerticesIndices.push_back(-1);
					}
				}
			}
		}
#ifndef SINGLE_THREADED_LOADING
#pragma omp single
#endif // !SINGLE_THREADED_LOADING
		{
			for(auto& face : surplusFaces) {
				bool finished = false;
				int startPos = (slotSize + 1) * face.first;
				int start = m_faceIndices[startPos];
				while(start != -1) {
					startPos = start;
					start = m_faceIndices[start];
				}
				bool added = false;
				for(int i = 0; i < slotSize; i++) {
					if(m_faceIndices[startPos + 1 + i] == -1) {
						m_faceIndices[startPos + 1 + i] = face.second;
						added = true;
						break;
					}
				}
				if(!added) {
					m_faceIndices[startPos] = static_cast<int>(m_faceIndices.size());
					m_faceIndices.push_back(-1);
					m_faceIndices.push_back(face.second);
					for(int i = 0; i < slotSize - 1; i++) {
						m_faceIndices.push_back(-1);
					}
				}
			}
		}
	}

	return true;
}

void ParallelMesh::loadFromData(std::vector<Vertex>&vertices, std::vector<int>&adjacentVerticesIndices, std::vector<int>&faceIndices, std::vector<glm::ivec3>&faces)
{
	m_faces = faces;
	m_vertices = vertices;
	m_adjacentVerticesIndices = adjacentVerticesIndices;
	m_faceIndices = faceIndices;
#ifdef MULTI_QUEUE
	m_priorityStructure = std::make_unique<PRIORITY_STRUCTURE>(static_cast<int>(vertices.size()), NUM_THREADS * 2);
#else
	m_priorityStructure = std::make_unique<PRIORITY_STRUCTURE>(static_cast<int>(vertices.size()));
#endif
}

void ParallelMesh::exportData(std::vector<Vertex>&vertices, std::vector<int>&adjacentVerticesIndices, std::vector<int>&faceIndices, std::vector<glm::ivec3>&faces)
{
	faces = m_faces;
	vertices = m_vertices;
	adjacentVerticesIndices = m_adjacentVerticesIndices;
	faceIndices = m_faceIndices;
}

std::pair<std::vector<glm::vec3>, std::vector<glm::vec3>> ParallelMesh::getData() {
	std::vector<glm::vec3> normals(m_vertices.size());
	std::vector<int> normalCount(m_vertices.size());

	std::vector<glm::ivec3> faces = this->getFaces();

	for (auto& face : faces) {
		glm::vec3 faceNormal = glm::cross(m_vertices[face.y].getPos() - m_vertices[face.x].getPos(), m_vertices[face.z].getPos() - m_vertices[face.x].getPos()); // right direction?
		normals[face.x] += faceNormal;
		normalCount[face.x] += 1;
		normals[face.y] += faceNormal;
		normalCount[face.y] += 1;
		normals[face.z] += faceNormal;
		normalCount[face.z] += 1;
	}
	for (int i = 0; i < normals.size(); i++)
	{

		normals[i] /= normalCount[i];
		normals[i] = glm::normalize(normals[i]);
		if (normals[i] != normals[i])
			normals[i] = { 1,0,0 };
	}

	std::vector<glm::vec3> normalsData;
	std::vector<glm::vec3> vertexData;
	for (auto& face : faces) {
		normalsData.push_back(normals[face.x]);
		normalsData.push_back(normals[face.y]);
		normalsData.push_back(normals[face.z]);
		vertexData.push_back(m_vertices[face.x].getPos());
		vertexData.push_back(m_vertices[face.y].getPos());
		vertexData.push_back(m_vertices[face.z].getPos());
	}


	return std::pair<std::vector<glm::vec3>, std::vector<glm::vec3>>(vertexData, normalsData);
}
std::vector<glm::vec3> ParallelMesh::getVertices()
{
	std::vector<glm::vec3> vertices(m_vertices.size());
	for (int i = 0; i < m_vertices.size(); i++) 
	{
		vertices[i] = m_vertices[i].getPos();
	}
	return vertices;
}
std::vector<glm::ivec3> ParallelMesh::getFaces() {
	std::vector<glm::ivec3> retFaces;
	for (int i = 0; i < m_vertices.size(); i++)
	{
		return m_faces;
	}
	return retFaces;
}

void ParallelMesh::computeQuadricErrorMatrices(int countThreads) {
	std::vector<std::pair<float, int>> errorValues(m_vertices.size());
#pragma omp parallel for num_threads(countThreads) schedule(guided)
	for (int i = 0; i < m_vertices.size(); i++)
	{	
		m_vertices[i].setErrorMetric();
		errorValues[i] = (std::pair<float, int>(m_vertices[i].getErrorMetric().getSmallestError(), i));
	}
#ifndef SINGLE_THREADED
	m_priorityStructure->setErrors(errorValues, countThreads);
#else
	m_priorityStructure->setErrors(errorValues);
#endif // !SINGLE_THREADED

}

float ParallelMesh::getVertexError(int index)
{
	return m_vertices[index].getErrorMetric().getSmallestError();
}

void ParallelMesh::updatePriorityStructure(int vertexId)
{
	m_priorityStructure->update(vertexId);
}

void ParallelMesh::debugCheckData(int activeException)
{
	int num_threads = 32;
//#pragma omp parallel for num_threads(num_threads)
	for (int i = 0; i < m_vertices.size(); i++) {
		if (m_vertices[i].isActive()) {
			std::vector<int> foundVertices = std::vector<int>();
			std::vector<bool> isValidAdjacentVertex = std::vector<bool>();
			{ 	//remove/replace Vertex from adjacentVertexList
				auto& adjacentVertices = ParallelMesh::getInstance().m_adjacentVerticesIndices;
				int startPos = (slotSize + 1) * i;
				int next;// = adjacentVertices[startPos];
				bool finished = false;
				while (startPos != -1) {
					next = adjacentVertices[startPos];
					for (int j = 0; j < slotSize; j++) {
						if (!finished) {
							if (std::find(foundVertices.begin(), foundVertices.end(), adjacentVertices[startPos + 1 + j]) != foundVertices.end()) {
								__debugbreak(); // Element multiple times in vector
							}
							else {
								if (adjacentVertices[startPos + 1 + j] != -1) {
									if (activeException != adjacentVertices[startPos + 1 + j] && !m_vertices[adjacentVertices[startPos + 1 + j]].isActive())// nonActive Vertex in list
										__debugbreak();
									foundVertices.push_back(adjacentVertices[startPos + 1 + j]);
									isValidAdjacentVertex.push_back(false);
								}
								else
									finished = true;
							}
						}
						else {
							if (adjacentVertices[startPos + 1 + j] != -1)
								__debugbreak(); // not clean
						}
					}
					startPos = next;
				}
			}			

			{
				auto& faceIndices = ParallelMesh::getInstance().m_faceIndices;
				int startPos = (3*slotSize + 1) * i;
				/*if (startPos == 4654)
					__debugbreak();*/
				int next;
				bool finished = false;
				int numFaces = 0;
				while (startPos != -1) {
					next = faceIndices[startPos];
					for (int j = 0; j < slotSize; j++) {
						if (!finished) {							
							if (faceIndices[startPos + 1 + 3 * j] != -1) {
								if (faceIndices[startPos + 1 + 3 * j + 1] == -1)
									__debugbreak(); // not clean
								else if (faceIndices[startPos + 1 + 3 * j + 2] == -1)
									__debugbreak(); // not clean
								else if (faceIndices[startPos + 1 + 3 * j] == faceIndices[startPos + 1 + 3 * j + 1] || faceIndices[startPos + 1 + 3 * j] == faceIndices[startPos + 1 + 3 * j + 2] || faceIndices[startPos + 1 + 3 * j + 2] == faceIndices[startPos + 1 + 3 * j + 1])
									__debugbreak();
								numFaces++;
								auto foundVertex = std::find(foundVertices.begin(), foundVertices.end(), faceIndices[startPos + 1 + 3 * j]);
								if (foundVertex != foundVertices.end()) {
									isValidAdjacentVertex[foundVertex - foundVertices.begin()] = true;
								}
								else {
									__debugbreak();
								}
								foundVertex = std::find(foundVertices.begin(), foundVertices.end(), faceIndices[startPos + 1 + 3 * j + 1]);
								if (foundVertex != foundVertices.end()) {
									isValidAdjacentVertex[foundVertex - foundVertices.begin()] = true;
								}
								else {
									__debugbreak();
								}
								foundVertex = std::find(foundVertices.begin(), foundVertices.end(), faceIndices[startPos + 1 + 3 * j + 2]);
								if (foundVertex != foundVertices.end()) {
									isValidAdjacentVertex[foundVertex - foundVertices.begin()] = true;
								}
								else {
									__debugbreak();
								}
							}
							else
								finished = true;
							
						}
						else {
							if (faceIndices[startPos + 1 + 3 * j] != -1)
								__debugbreak(); // not clean
							if (faceIndices[startPos + 1 + 3 * j + 1] != -1)
								__debugbreak(); // not clean
							if (faceIndices[startPos + 1 + 3 * j + 2] != -1)
								__debugbreak(); // not clean
						}
					}
					startPos = next;
				}
				if (numFaces < 3)
					__debugbreak();
			}
			for (bool b : isValidAdjacentVertex) {
				if (!b)
					__debugbreak();
			}
		}
	}
}

void ParallelMesh::addTriangle(int currentIndex, int index1, int index2, int a, int b, int c, CRITICAL_SECTION& csF, CRITICAL_SECTION& csV, std::vector<std::pair<int, int>>& surplusVertices, std::vector<std::pair<int, glm::ivec3>>& surplusFaces)
{
	bool found1, found2, foundCurrent;
	int i;
	if (m_faceIndices[(currentIndex + 1) * (3 * slotSize + 1) - 1] != -1) {
		EnterCriticalSection(&csF);
		surplusFaces.push_back({ currentIndex,{a,b,c} });
		LeaveCriticalSection(&csF);
	}
	else {
		for (int i = 0; i < slotSize; i++) {
			if (m_faceIndices[currentIndex * (3 * slotSize + 1) + 1 + 3 * i] == -1) {
				m_faceIndices[currentIndex * (3 * slotSize + 1) + 1 + 3 * i] = a;
				m_faceIndices[currentIndex * (3 * slotSize + 1) + 1 + 3 * i + 1] = b;
				m_faceIndices[currentIndex * (3 * slotSize + 1) + 1 + 3 * i + 2] = c;
				break;
			}
		}
	}
	foundCurrent = false;
	found1 = false;
	found2 = false;
	for (int i = 0; i < slotSize; i++) {
		if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == currentIndex) {
			foundCurrent = true;
		}
		if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == index1) {
			found1 = true;
		}
		else if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == index2) {
			found2 = true;
		}
	}
	i = 0;
	if (!foundCurrent) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex, currentIndex });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = currentIndex;
					break;
				}
			}
		}
	}
	if (!found1) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex,index1 });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = index1;
					break;
				}
			}
		}		
	}
	if (!found2) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) -1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex,index2 });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = index2;
					break;
				}
			}
		}		
	}
}

void ParallelMesh::addTriangleIndex(int triangleIndex, int currentIndex, int index1, int index2, int a, int b, int c, CRITICAL_SECTION& csF, CRITICAL_SECTION& csV, std::vector<std::pair<int, int>>& surplusVertices, std::vector<std::pair<int, int>>& surplusFaces)
{
	bool found1, found2, foundCurrent;
	if (m_faceIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
		EnterCriticalSection(&csF);
		surplusFaces.push_back({ currentIndex,triangleIndex });
		LeaveCriticalSection(&csF);
	}
	else {
		for (int i = 0; i < slotSize; i++) {
			if (m_faceIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
				m_faceIndices[currentIndex * (slotSize + 1) + 1 + i] = triangleIndex;
				break;
			}
		}
	}
	foundCurrent = false;
	found1 = false;
	found2 = false;
	for (int i = 0; i < slotSize; i++) {
		if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == currentIndex) {
			foundCurrent = true;
		}
		if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == index1) {
			found1 = true;
		}
		else if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == index2) {
			found2 = true;
		}
	}
	int i = 0;
	if (!foundCurrent) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex, currentIndex });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = currentIndex;
					break;
				}
			}
		}
	}
	if (!found1) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex,index1 });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = index1;
					break;
				}
			}
		}
	}
	if (!found2) {
		if (m_adjacentVerticesIndices[(currentIndex + 1) * (slotSize + 1) - 1] != -1) {
			EnterCriticalSection(&csV);
			surplusVertices.push_back({ currentIndex,index2 });
			LeaveCriticalSection(&csV);
		}
		else {
			for (; i < slotSize; i++) {
				if (m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] == -1) {
					m_adjacentVerticesIndices[currentIndex * (slotSize + 1) + 1 + i] = index2;
					break;
				}
			}
		}
	}
}


void ParallelMesh::reduceVertices(int finalSize, int countThreads) {

	std::atomic<int> indexCounter = 0;
#ifdef VECTOR_REDUCTION
	std::vector<int> vectorInit(m_vertices.size());
#ifndef SINGLE_THREADED
	initReplacementStart(countThreads);
#pragma omp parallel num_threads(countThreads)
	{
		int threadNr = omp_get_thread_num();
		setReplacementStart(m_priorityStructure->getSize(), countThreads, threadNr);
#pragma omp barrier
		//std::sort(std::execution::par_unseq, nodes.begin(), nodes.begin() + finalSize - 1, std::less<int>());

		compact(m_priorityStructure->getSize(), finalSize, countThreads, threadNr, vectorInit, m_priorityStructure->getNodes());
#pragma omp barrier
	}
#else
	compact(m_vertices.size(), finalSize, 1, 0, vectorInit, m_priorityStructure->getNodes());
#endif

//#endif // !SINGLE_THREADED
	const std::vector<int> map(std::move(vectorInit));
#else
	std::unordered_map<int, int> map(finalSize);
	for (int i = 0; i < m_vertices.size(); i++)
	{
		if (m_vertices[i].isActive())
		{
			map.insert(std::pair<int, int>(i, indexCounter));
			m_vertices[indexCounter++] = m_vertices[i];
		}
	}
#endif // !VECTOR_REDUCTION
	m_vertices.erase(m_vertices.begin() + finalSize, m_vertices.end());

	auto& faceIndices = m_faceIndices;
	auto& adjacentVerticesIndices = m_adjacentVerticesIndices;

	std::vector<int> newAdjacentVerticesIndices;
	std::vector<int> newFaceIndices;

	auto& faces = m_faces;

	std::vector<glm::ivec3> newFaces(std::min(int(m_vertices.size() * 2.1) + 1000, int(m_faces.size())));
	std::vector<int> facesMap(m_faces.size(), -1);
#ifndef SINGLE_THREADED
	std::atomic<int> FOverflow, VOverflow;
	std::vector<int> FOverflowNr(m_vertices.size());
	std::vector<int> VOverflowNr(m_vertices.size());
	std::vector<std::atomic<int>> facesMapLocks(m_faces.size());
#endif // !SINGLE_THREADED

	std::atomic<int> faceCounter = 0;
#ifndef SINGLE_THREADED
#pragma omp parallel num_threads(countThreads)
	{
#pragma omp for schedule(guided)
		for (int i = 0; i < m_vertices.size(); i++) {
			{
				int overflowCount = 0;
				int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
				int readPos = readStart + 1;
				bool stop = false;

				while (!stop) {
					if (faceIndices[readStart] != -1)
					{
						overflowCount++;
						readStart = faceIndices[readStart];
					}
					else {
						stop = true;
					}
				}
				if (overflowCount) {
					FOverflowNr[i] = (FOverflow += overflowCount);
					FOverflowNr[i] -= overflowCount;
				}
			}
			{
				int overflowCount = 0;
				int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
				int readPos = readStart + 1;
				bool stop = false;

				while (!stop) {
					if (adjacentVerticesIndices[readStart] != -1)
					{
						overflowCount++;
						readStart = adjacentVerticesIndices[readStart];
					}
					else {
						stop = true;
					}
				}
				if (overflowCount) {
					VOverflowNr[i] = (VOverflow += overflowCount);
					VOverflowNr[i] -= overflowCount;
				}
			}
		}
		if (omp_get_thread_num() == 0)
		{
			int sizeVertices = (slotSize + 1) * (static_cast<int>(m_vertices.size()) + VOverflow);
			int sizeFaces = (slotSize + 1) * (static_cast<int>(m_vertices.size()) + FOverflow);
			newAdjacentVerticesIndices = std::vector<int>(sizeVertices, -1);
			newFaceIndices = std::vector<int>(sizeFaces, -1);
		}
#pragma omp barrier
#pragma omp for schedule(guided)
		for (int i = 0; i < m_vertices.size(); i++)
		{
			//if (m_vertices[i].isActive())
			{
				int overflowCount = 0;
				int writeStart = (slotSize + 1) * i;
				int writePos = writeStart + 1;
				int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
				int readPos = readStart + 1;
				bool stop = false;
				while (faceIndices[readPos] != -1 && !stop)
				{
					for (int j = 0; j < slotSize; j++) {
						if (faceIndices[readPos] != -1) {
#ifndef VECTOR_REDUCTION
							__debugbreak();
#else
							int faceIndex = facesMap[faceIndices[readPos]];
							if (faceIndex == -1)
							{
#ifndef SINGLE_THREADED
								if (facesMapLocks[faceIndices[readPos]]++ == 0)
								{
#endif // !SINGLE_THREADED
									if (facesMap[faceIndices[readPos]] == -1) {
										faceIndex = faceCounter++;
										newFaces[faceIndex] = { map[faces[faceIndices[readPos]].x], map[faces[faceIndices[readPos]].y], map[faces[faceIndices[readPos]].z] };
										facesMap[faceIndices[readPos]] = faceIndex;
									}
									else {
										faceIndex = facesMap[faceIndices[readPos]];
									}
#ifndef SINGLE_THREADED
									facesMapLocks[faceIndices[readPos]]--;
								}
								else {
									facesMapLocks[faceIndices[readPos]]--;
									while (facesMapLocks[faceIndices[readPos]] != 0);
									faceIndex = facesMap[faceIndices[readPos]];
								}
#endif
							}
							readPos++;
							newFaceIndices[writePos++] = faceIndex;
#endif // VECTOR_REDUCTION
						}// else break?
					}
					if (faceIndices[readStart] != -1)
					{
						readStart = faceIndices[readStart];
						readPos = readStart + 1;
						int newWriteStart = (static_cast<int>(m_vertices.size()) + FOverflowNr[i] + overflowCount) * (slotSize + 1);
						newFaceIndices[writeStart] = newWriteStart;
						writeStart = newWriteStart;
						writePos = writeStart + 1;
						overflowCount++;
					}
					else {
						stop = true;
					}

				}
			}
			{
				int overflowCount = 0;
				int writeStart = (slotSize + 1) * i;
				int writePos = writeStart + 1;
				int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
				int readPos = readStart + 1;
				bool stop = false;

				while (adjacentVerticesIndices[readPos] != -1 && !stop)
				{
					for (int j = 0; j < slotSize; j++) {
						if (adjacentVerticesIndices[readPos] != -1) {
#ifndef VECTOR_REDUCTION
							newAdjacentVerticesIndices[writePos++] = map[adjacentVerticesIndices[readPos++]]; // Copy everything from this slot
#else
							newAdjacentVerticesIndices[writePos++] = map[adjacentVerticesIndices[readPos++]];
#endif // VECTOR_REDUCTION
						}
						if (adjacentVerticesIndices[readStart] != -1)
						{
							readStart = adjacentVerticesIndices[readStart];
							readPos = readStart + 1;
							int newWriteStart = (static_cast<int>(m_vertices.size()) + VOverflowNr[i] + overflowCount) * (slotSize + 1);
							newAdjacentVerticesIndices[writeStart] = newWriteStart;
							writeStart = newWriteStart;
							writePos = writeStart + 1;
							overflowCount++;
						}
						else {
							stop = true;
						}
					}
				}
			}
		}
#pragma omp barrier
		if (omp_get_thread_num() == 0)
		{
			newFaces.erase(newFaces.begin() + faceCounter, newFaces.end());
		}
	}
#else
	int VOverflow = 0;
	int FOverflow = 0;
	for (int i = 0; i < m_vertices.size(); i++) {
		{
			int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
			int readPos = readStart + 1;
			bool stop = false;

			while (!stop) {
				if (faceIndices[readStart] != -1)
				{
					FOverflow++;
					readStart = faceIndices[readStart];
				}
				else {
					stop = true;
				}
			}
		}
		{
			int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
			int readPos = readStart + 1;
			bool stop = false;

			while (!stop) {
				if (adjacentVerticesIndices[readStart] != -1)
				{
					VOverflow++;
					readStart = adjacentVerticesIndices[readStart];
				}
				else {
					stop = true;
				}
			}

		}
	}

	int sizeVertices = (slotSize + 1) * (m_vertices.size() + VOverflow);
	int sizeFaces = (slotSize + 1) * (m_vertices.size() + FOverflow);
	newAdjacentVerticesIndices = std::vector<int>(sizeVertices, -1);
	newFaceIndices = std::vector<int>(sizeFaces, -1);

	int overflowCountF = 0;
	int overflowCountV = 0;
	for (int i = 0; i < m_vertices.size(); i++)
	{
		//if (m_vertices[i].isActive())
		{
			int writeStart = (slotSize + 1) * i;
			int writePos = writeStart + 1;
			int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
			int readPos = readStart + 1;
			bool stop = false;
			while (faceIndices[readPos] != -1 && !stop)
			{
				for (int j = 0; j < slotSize; j++) {
					if (faceIndices[readPos] != -1) {
#ifndef VECTOR_REDUCTION
						__debugbreak();
#else
						int faceIndex = facesMap[faceIndices[readPos]];
						if (faceIndex == -1)
						{
							if (facesMap[faceIndices[readPos]] == -1) {
								faceIndex = faceCounter++;
								newFaces[faceIndex] = { map[faces[faceIndices[readPos]].x], map[faces[faceIndices[readPos]].y], map[faces[faceIndices[readPos]].z] };
								facesMap[faceIndices[readPos]] = faceIndex;
							}
							else {
								faceIndex = facesMap[faceIndices[readPos]];
							}
#ifndef SINGLE_THREADED
#endif
						}
						readPos++;
						newFaceIndices[writePos++] = faceIndex;
#endif // VECTOR_REDUCTION
					}// else break?
				}
				if (faceIndices[readStart] != -1)
				{
					readStart = faceIndices[readStart];
					readPos = readStart + 1;
					int newWriteStart = (m_vertices.size() + overflowCountF) * (slotSize + 1);
					newFaceIndices[writeStart] = newWriteStart;
					writeStart = newWriteStart;
					writePos = writeStart + 1;
					overflowCountF++;
				}
				else {
					stop = true;
				}

			}
		}
		{
			int writeStart = (slotSize + 1) * i;
			int writePos = writeStart + 1;
			int readStart = (slotSize + 1) * m_vertices[i].getVertexIndex();
			int readPos = readStart + 1;
			bool stop = false;

			while (adjacentVerticesIndices[readPos] != -1 && !stop)
			{
				for (int j = 0; j < slotSize; j++) {
					if (adjacentVerticesIndices[readPos] != -1) {
#ifndef VECTOR_REDUCTION
						newAdjacentVerticesIndices[writePos++] = map[adjacentVerticesIndices[readPos++]]; // Copy everything from this slot
#else
						newAdjacentVerticesIndices[writePos++] = map[adjacentVerticesIndices[readPos++]];
#endif // VECTOR_REDUCTION
					}
					if (adjacentVerticesIndices[readStart] != -1)
					{
						readStart = adjacentVerticesIndices[readStart];
						readPos = readStart + 1;
						int newWriteStart = (m_vertices.size() + overflowCountV) * (slotSize + 1);
						newAdjacentVerticesIndices[writeStart] = newWriteStart;
						writeStart = newWriteStart;
						writePos = writeStart + 1;
						overflowCountV++;
					}
					else {
						stop = true;
					}
				}
			}
		}

	}
	newFaces.erase(newFaces.begin() + faceCounter, newFaces.end());

#endif

	m_faces.swap(newFaces);
	m_adjacentVerticesIndices.swap(newAdjacentVerticesIndices);
	m_faceIndices.swap(newFaceIndices);
}
