#include "Vertex.h"
#include "ParallelMesh.h"

Vertex::Vertex(glm::vec3 pos, int index) : m_pos(pos), m_active(true), m_vertexIndex(index)
{	
	//m_lock.clear();
}

void Vertex::changePos(glm::vec3 pos) {
	m_pos = pos;
}
bool Vertex::deactivate() {
	// returns true if is not already deactivated;
	bool returnValue = m_active;
	m_active = false;
	return returnValue;
}
bool Vertex::isActive() {
	return m_active;
}


bool Vertex::removeTriangle(int triangleIndex) {
	auto& faceIndices = ParallelMesh::getInstance().m_faceIndices;

	int startPos = (slotSize + 1) * m_vertexIndex;
	int lastStartPos = startPos;
	int start = faceIndices[startPos];
	int currentPos = startPos + 1;
	int ctr = 0;
	bool removed = false;
	bool stop = false;
	while (!stop && faceIndices[currentPos] != -1) {
		if (faceIndices[currentPos] == triangleIndex) {
			// find last face to replace this face
			if (start == -1 && (ctr == slotSize - 1 || faceIndices[currentPos + 1] == -1)) { // is this the last face?
				if (ctr == 0) {
					faceIndices[currentPos] = -1;
					faceIndices[lastStartPos] = -1;
				}
				else {
					faceIndices[currentPos] = -1;
					faceIndices[startPos] = -1;
				}

			}
			else {
				int lastIndex = currentPos;

				while (start != -1) {
					lastStartPos = startPos;
					startPos = start;
					start = faceIndices[startPos];
					lastIndex = startPos + 1;
					ctr = 0;
				}

				while (++ctr != slotSize && faceIndices[lastIndex + 1] != -1) {
					//++ctr;
					lastIndex++;
				}
				faceIndices[currentPos] = faceIndices[lastIndex];
				faceIndices[lastIndex] = -1;
				if (ctr == 1) {
					faceIndices[lastStartPos] = -1;
				}
			}
			/*if (ctr == 0) {
				faceIndices[startPos] = -1;
			}*/
			removed = true;
			break;
		}
		if (++ctr != slotSize) {
			currentPos++;
		}
		else {
			ctr = 0;
			if (start != -1)
			{
				lastStartPos = startPos;
				startPos = start;
				start = faceIndices[startPos];
				currentPos = startPos + 1;
			}
			else {
				stop = true;
			}
		}
	}
	if (removed) {
		volatile int pos = (slotSize + 1) * m_vertexIndex + 1;
		if (faceIndices[pos] == -1) { // Vertex with no Triangles
			//__debugbreak();
			//m_adjacentVerticesIndices.clear();
			//updateErrorMetric();
			//return true;
		}
		return false;
	}

	throw new std::invalid_argument("Non existing triangle deleted");
}

void Vertex::setErrorMetric() {
	m_error = QuadricErrorMetric();
	m_error.init(*this);
}
QuadricErrorMetric& Vertex::getErrorMetric()
{
	return m_error;
}
void Vertex::updateErrorMetric(int newVertex)
{
	float oldError = m_error.getSmallestError();
	if(newVertex == -1)
		m_error.update(*this);
	else
		m_error.updateSingle(newVertex);
	float newError = m_error.getSmallestError();
	if (newError != oldError) 
	{
		ParallelMesh::getInstance().updatePriorityStructure(m_vertexIndex);
	}
}

void Vertex::collapse(int ownIndex)
{
#ifndef NO_DEBUG_CHECKS
	if (!this->deactivate()) {
		__debugbreak();
	}
#else
	this->deactivate();
#endif // !NO_DEBUG_CHECKS
	auto& parent = ParallelMesh::getInstance();
	int vertexIndexToCollapseTo = m_error.getSmallestErrorIndex();
	if (vertexIndexToCollapseTo == ownIndex) {
#ifndef SINGLE_THREADED
		parent.decreaseAdjacentCollapses(m_vertexIndex);
#endif
		return;
	}
#ifndef NO_DEBUG_CHECKS
	if (vertexIndexToCollapseTo < 0) {
		// We have no faces
		__debugbreak();
	}	
	else 
#endif // !NO_DEBUG_CHECKS
	{		

		auto& adjacentVertices = parent.m_adjacentVerticesIndices;

		{
			int startPos = (slotSize + 1) * m_vertexIndex;
			int start = adjacentVertices[startPos];
			int currentPos = startPos + 1;
			int ctr = 0;
			bool stop = false;
			parent.getVertex(vertexIndexToCollapseTo).collapseUpdateTriangles(ownIndex, vertexIndexToCollapseTo); // It is important to do this first so the changed triangles are valid for all surrounding vertices when they are released
			while (!stop && adjacentVertices[currentPos] != -1) {

				if ((adjacentVertices[currentPos] != ownIndex) && (adjacentVertices[currentPos] != vertexIndexToCollapseTo)) {
					parent.getVertex(adjacentVertices[currentPos]).incommingCollapse(ownIndex, vertexIndexToCollapseTo, false); // Replace Vertex
#ifndef SINGLE_THREADED
					parent.decreaseAdjacentCollapses(adjacentVertices[currentPos]); // We dont work with this vertex Anymore so we can release it early
#endif
				}

				if (++ctr != slotSize) {
					currentPos += 1;
				}
				else {
					ctr = 0;
					if (start != -1)
					{
						startPos = start;
						start = adjacentVertices[startPos];
						currentPos = startPos + 1;
					}
					else {
						stop = true;
					}
				}
			}
		}
		parent.getVertex(vertexIndexToCollapseTo).incommingCollapse(ownIndex, vertexIndexToCollapseTo,true); // Its important to do this last because otherwise the adjacentVertices are not valid anymore
		
#ifndef SINGLE_THREADED
		parent.decreaseAdjacentCollapses(vertexIndexToCollapseTo); // We dont work with this vertex Anymore so we can release it.
		//parent.decreaseAdjacentCollapses(m_vertexIndex); // Vertex is Deleted
#endif

	}
}

bool Vertex::isValidCollapse()
{
	// Check if the vertices have exactly two common neigbours
	auto& adjacentVertices = ParallelMesh::getInstance().m_adjacentVerticesIndices;
	int countCommonNeigbours = 0;
	int collapseVertexIndex = m_error.getSmallestErrorIndex();
#ifndef SINGLE_THREADED
	//ParallelMesh::getInstance().lockVertex(collapseVertexIndex);
#endif
	{
		int startPos = (slotSize + 1) * m_vertexIndex;
		int start = adjacentVertices[startPos];
		int currentPos1 = startPos + 1;
		int ctr = 0;
		bool stop = false;
		while (!stop && adjacentVertices[currentPos1] != -1) {

			{
				int startPos = (slotSize + 1) * collapseVertexIndex;
				int start = adjacentVertices[startPos];
				int currentPos = startPos + 1;
				int ctr = 0;
				bool stop = false;
				while (!stop && adjacentVertices[currentPos] != -1) {

					if (adjacentVertices[currentPos] == adjacentVertices[currentPos1])
					{
						countCommonNeigbours++;
						break;
					}


					if (++ctr != slotSize) {
						currentPos += 1;
					}
					else {
						ctr = 0;
						if (start != -1)
						{
							startPos = start;
							start = adjacentVertices[startPos];
							currentPos = startPos + 1;
						}
						else {
							stop = true;
						}
					}
				}
			}

			
			if (++ctr != slotSize) {
				currentPos1 += 1;
			}
			else {
				ctr = 0;
				if (start != -1)
				{
					startPos = start;
					start = adjacentVertices[startPos];
					currentPos1 = startPos + 1;
				}
				else {
					stop = true;
				}
			}
		}
	}
#ifndef SINGLE_THREADED
	//ParallelMesh::getInstance().unlockVertex(collapseVertexIndex);
#endif
	if (countCommonNeigbours != 4) { // itselfs + 2 neigbours are allowed
		//__debugbreak();
		return false;
	}

#ifdef NORMAL_VALIDATION_CHECK
	auto& parent = ParallelMesh::getInstance();
	auto& faceIndices = parent.m_faceIndices;
	auto& faces = parent.m_faces;
	int startPos = (slotSize + 1) * m_vertexIndex;
	int currentPos = startPos + 1;
	int ctr = 0;
	bool stop = false;
	glm::vec3 otherNormal = parent.getVertex(collapseVertexIndex).getPos();
	while (!stop && faceIndices[currentPos] != -1) {
		{
			if (faces[faceIndices[currentPos]].x == m_vertexIndex) {
				if (faces[faceIndices[currentPos]].y == collapseVertexIndex || faces[faceIndices[currentPos]].z == collapseVertexIndex) {

				}
				else {
					glm::vec3 v0, v1, v2;
					v0 = m_pos;
					v1 = parent.getVertex(faces[faceIndices[currentPos]].y).getPos();
					v2 = parent.getVertex(faces[faceIndices[currentPos]].z).getPos();
					glm::dvec3 faceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					v0 = otherNormal;
					glm::dvec3 otherFaceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					if (glm::dot(otherFaceNormal, faceNormal) < NORMAL_THRESHOLD)
						return false;
				}
			}
			else if (faces[faceIndices[currentPos]].y == m_vertexIndex) {
				if (faces[faceIndices[currentPos]].x == collapseVertexIndex || faces[faceIndices[currentPos]].z == collapseVertexIndex) {

				}
				else {
					glm::vec3 v0, v1, v2;
					v0 = m_pos;
					v1 = parent.getVertex(faces[faceIndices[currentPos]].x).getPos();
					v2 = parent.getVertex(faces[faceIndices[currentPos]].z).getPos();
					glm::dvec3 faceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					v0 = otherNormal;
					glm::dvec3 otherFaceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					if (glm::dot(otherFaceNormal, faceNormal) < NORMAL_THRESHOLD)
						return false;
				}
			}
			else if (faces[faceIndices[currentPos]].z == m_vertexIndex) {
				if (faces[faceIndices[currentPos]].x == collapseVertexIndex || faces[faceIndices[currentPos]].y == collapseVertexIndex) {

				}
				else {
					glm::vec3 v0, v1, v2;
					v0 = m_pos;
					v1 = parent.getVertex(faces[faceIndices[currentPos]].y).getPos();
					v2 = parent.getVertex(faces[faceIndices[currentPos]].x).getPos();
					glm::dvec3 faceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					v0 = otherNormal;
					glm::dvec3 otherFaceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
					if (glm::dot(otherFaceNormal, faceNormal) < NORMAL_THRESHOLD)
						return false;
				}
			}
		}
		if (++ctr != slotSize) {
			currentPos++;
		}
		else {
			if (faceIndices[startPos] != -1)
			{
				ctr = 0;
				startPos = faceIndices[startPos];
				//start = faceIndices[startPos];
				currentPos = startPos + 1;
			}
			else {
				ctr = slotSize - 1;
				stop = true;
			}
		}
	}
#endif
	return true;
}

void Vertex::incommingCollapse(int indexCollapse, int indexToCollapseTo, bool isIndexToCollapseTo)
{
#ifndef NO_DEBUG_CHECKS
	if (!this->isActive())
		__debugbreak();
#endif // !NO_DEBUG_CHECKS
	auto& adjacentVertices = ParallelMesh::getInstance().m_adjacentVerticesIndices;
	int replacementPos = -1;
	{
		while (replacementPos == -1)
		{ 	//remove/replace Vertex from adjacentVertexList
			int startPos = (slotSize + 1) * m_vertexIndex;
			int start = adjacentVertices[startPos];
			int currentPos = startPos + 1;
			int ctr = 0;
			bool stop = false;

			int lastPos = -1;
			bool indexToCollapseToFound = false;
			while (!stop && adjacentVertices[currentPos] != -1) {
				if (adjacentVertices[currentPos] == indexCollapse)
				{
					replacementPos = currentPos;
				}
				else if (adjacentVertices[currentPos] == indexToCollapseTo) {
					indexToCollapseToFound = true;
				}
				lastPos = currentPos;
				if (++ctr != slotSize) {
					currentPos += 1;
				}
				else {
					ctr = 0;
					if (start != -1)
					{
						startPos = start;
						start = adjacentVertices[startPos];
						currentPos = startPos + 1;
					}
					else {
						stop = true;
					}
				}
			}
#ifndef NO_DEBUG_CHECKS
			if (replacementPos == -1) {
				__debugbreak();		
			}
			else
#endif // !NO_DEBUG_CHECKS
			{
				if (!indexToCollapseToFound) {
					adjacentVertices[replacementPos] = indexToCollapseTo;
				}
				else { // swap with last element
					adjacentVertices[replacementPos] = adjacentVertices[lastPos];
					adjacentVertices[lastPos] = -1; // also works when lastPos == replacementPos
				}
			}
		}
	}
	if(isIndexToCollapseTo) {
		// Merge adjacent Vertices

		int startPosStealV = (slotSize + 1) * indexCollapse;;
		int stolenStartPos;
		bool stolenV = false;
		int lastLegitPosV = -1;

		{
			int startPos = (slotSize + 1) * indexCollapse;
			int start = adjacentVertices[startPos];
			int currentPos = startPos + 1;
			int ctr = 0;
			bool stop = false;

			while((!stop && adjacentVertices[currentPos] != -1 && currentPos != lastLegitPosV)) {
				if(adjacentVertices[currentPos] != indexCollapse) {
					int vertexToFind = adjacentVertices[currentPos];
					int startPos = (slotSize + 1) * m_vertexIndex;
					int start = adjacentVertices[startPos];
					int currentPos = startPos + 1;
					int ctr = 0;
					bool stop = false;
					bool found = false;
					bool legit = true;
					while(!stop && adjacentVertices[currentPos] != -1 && legit && !found) {
						if(currentPos == lastLegitPosV) {
							legit = false;
						}
						if(adjacentVertices[currentPos] == vertexToFind) {
							found = true;
						}
						if(++ctr != slotSize) {
							currentPos += 1;
						} else {
							ctr = 0;
							if(start != -1) {
								startPos = start;
								start = adjacentVertices[startPos];
								currentPos = startPos + 1;
							} else {
								stop = true;
							}
						}
					}
					if(!found) {
						if(stop) {
							// steal the storage from the vertex to delete
							adjacentVertices[startPos] = startPosStealV;
							stolenStartPos = startPosStealV;
							currentPos = startPosStealV + 1;
							startPosStealV = adjacentVertices[startPosStealV];
							stolenV = true;
						} else {
#ifndef NO_DEBUG_CHECKS
							if(adjacentVertices[currentPos] != -1 && !stolenV)
								__debugbreak();
#endif
						}
						adjacentVertices[currentPos] = vertexToFind;
						if(stolenV) {
							lastLegitPosV = currentPos;
						}
					}
				}
				if(++ctr != slotSize) {
					currentPos += 1;
				} else {
					ctr = 0;
					if(start != -1) {
						startPos = start;
						start = adjacentVertices[startPos];
						currentPos = startPos + 1;
					} else {
						stop = true;
					}
				}
			}
		}
		if(stolenV) {
			{
				int lastLegitCtr = lastLegitPosV % (slotSize + 1);
				int lastLegitStartPos = lastLegitPosV - lastLegitCtr;
				adjacentVertices[lastLegitStartPos] = -1;

				for(int i = 1; i - 1 < slotSize - lastLegitCtr; i++) {
					adjacentVertices[lastLegitPosV + i] = -1;
				}
			}
		}
		this->getErrorMetric().updateQuadric(ParallelMesh::getInstance().getVertex(indexCollapse).getErrorMetric().getQuadric());
		updateErrorMetric();
	}
#ifdef ERROR_UPDATE_OPTIMISED
	else if (this->getErrorMetric().getSmallestErrorIndex() == indexCollapse || this->getErrorMetric().getSmallestErrorIndex() == indexToCollapseTo) {
		updateErrorMetric();
	}
	else {
		updateErrorMetric(indexToCollapseTo);
	}
#else
	else {
		updateErrorMetric();
	}	
#endif
	if (NUM_THREADS == -1)
		ParallelMesh::getInstance().debugCheckData(indexCollapse);
}

void Vertex::collapseUpdateTriangles(int indexCollapse, int indexToCollapseTo) {
	auto& faceIndices = ParallelMesh::getInstance().m_faceIndices;
	auto& faces = ParallelMesh::getInstance().m_faces;
	int startPos = (slotSize + 1) * indexCollapse; // We could also iterate over our own faces
	int currentPos = startPos + 1;
	int ctr = 0;
	bool stop = false;
	while(!stop && faceIndices[currentPos] != -1) { // Remove the faces that are deleted in the HEC from the connected vertices (ignore vertex that is deleted)
		int triangle = faceIndices[currentPos];
		glm::ivec3 face = faces[faceIndices[currentPos]];
		// The faces that have to be deleted are connected to the collapsing vertex (no need to check we are iterating over its faces), indexToCollapseTo and an other vertex
		if(face.x == indexToCollapseTo || face.y == indexToCollapseTo || face.z == indexToCollapseTo) {
			// Remove the triangle from the other vertex;
			if(face.x != indexCollapse && face.x != indexToCollapseTo) {
				ParallelMesh::getInstance().getVertex(face.x).removeTriangle(triangle);
			} else if(face.y != indexCollapse && face.y != indexToCollapseTo) {
				ParallelMesh::getInstance().getVertex(face.y).removeTriangle(triangle);
			} else if(face.z != indexCollapse && face.z != indexToCollapseTo) {
				ParallelMesh::getInstance().getVertex(face.z).removeTriangle(triangle);
			}
			// Remove the triangle from ourself;
			removeTriangle(triangle);
		}
		if(++ctr != slotSize) {
			currentPos++;
		} else {
			if(faceIndices[startPos] != -1) {
				ctr = 0;
				startPos = faceIndices[startPos];
				currentPos = startPos + 1;
			} else {
				ctr = slotSize - 1;
				stop = true;
			}
		}
	}

	int startPosRead = (slotSize + 1) * indexCollapse;
	int startRead = faceIndices[startPosRead];
	int currentPosRead = startPosRead + 1;
	int ctrRead = 0;

	int startPosStealF = startPosRead;
	bool stolenF = false;

	int startPosWrite = (slotSize + 1) * m_vertexIndex;
	int ctrWrite = 0;
	bool stopWrite = false;
	while(faceIndices[startPosWrite] != -1) {
		startPosWrite = faceIndices[startPosWrite];
	}
	int currentPosWrite = startPosWrite + 1;
	while(!stopWrite && faceIndices[currentPosWrite] != -1) {
		currentPosWrite++;
		if(++ctrWrite == slotSize) {
			stopWrite = true;
		}
	}
	currentPosWrite--;
	ctrWrite--;

	// merge Faces
	stop = false;
	while(!stop && faceIndices[currentPosRead] != -1) {
		if(!(faces[faceIndices[currentPosRead]].x == indexToCollapseTo || faces[faceIndices[currentPosRead]].y == indexToCollapseTo || faces[faceIndices[currentPosRead]].z == indexToCollapseTo)) {
			bool matchesFirst = faces[faceIndices[currentPosRead]].x == indexCollapse;
			bool matchesSecond = faces[faceIndices[currentPosRead]].y == indexCollapse;
			bool matchesThird = faces[faceIndices[currentPosRead]].z == indexCollapse;
			if(matchesFirst || matchesSecond || matchesThird) {
				if(++ctrWrite != slotSize) { // We have still space in our current slot
					currentPosWrite++;
				} else {
					// steal the storage from the vertex to delete
					ctrWrite = 0;
					if(!stolenF)
						faceIndices[startPosWrite] = startPosStealF;
					stolenF = true;
#ifndef NO_DEBUG_CHECKS
					if(startPosStealF == -1)
						__debugbreak();
#endif // !NO_DEBUG_CHECKS
					currentPosWrite = startPosStealF + 1;
					startPosStealF = faceIndices[startPosStealF];
				}
				if(matchesFirst) {
					faces[faceIndices[currentPosRead]].x = indexToCollapseTo;
				} else if(matchesSecond) {
					faces[faceIndices[currentPosRead]].y = indexToCollapseTo;
				} else {
					faces[faceIndices[currentPosRead]].z = indexToCollapseTo;
				}
				faceIndices[currentPosWrite] = faceIndices[currentPosRead];

			}
#ifndef NO_DEBUG_CHECKS
			else {
				__debugbreak();
			}
#endif // !NO_DEBUG_CHECKS
		}

		if(++ctrRead != slotSize) {
			currentPosRead++;
		} else {
			ctrRead = 0;
			if(startRead != -1) {
				startPosRead = startRead;
				startRead = faceIndices[startPosRead];
				currentPosRead = startPosRead + 1;
			} else {
				stop = true;
			}
		}
	}

	if(stolenF) {
		while(++ctrWrite != slotSize) {
			currentPosWrite++;
			faceIndices[currentPosWrite] = -1;
		}
		faceIndices[currentPosWrite - (slotSize - 1) - 1] = -1;
	}
}
