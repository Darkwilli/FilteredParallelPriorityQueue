#pragma once

#include <glm/glm.hpp>
#include <vector>

#include "QuadricErrorMetric.h"

class ParallelMesh;
// Our Vertex class contains the Vertex Information: The position where it is located, a flag if it should be deleted because it was decimated, its own index and an QuadricErrorMetric
// The QEM contains the error quadric, the best vertex index to collapse to and its HEC error
class Vertex
{
public:
	friend ParallelMesh;
	friend ErrorMetric;
	friend QuadricErrorMetric;
	Vertex(){ };
	Vertex(glm::vec3 pos, int index);

	void changePos(glm::vec3 pos);
	glm::vec3 getPos() const;
	bool deactivate(); // Marks the vertex for deletion returns true if it is not already marked
	bool isActive();
	int getVertexIndex() const;
	bool removeTriangle(int triangleIndex); // removes the triangle from the vertex
	void setErrorMetric();
	QuadricErrorMetric& getErrorMetric();
	void updateErrorMetric(int newVertex = -1);
	float setNextBestError() { // calculates the next best error for the collapse (because the best is not a valid option) 
		return m_error.setNextBestError(*this);
	}
	// collapse the vertrex
	void collapse(int ownIndex);
	// Do One-Ring-Test and Normal test (if defined) to test if the current collapse (best according to error metric) is valid 
	bool isValidCollapse();

	// A neighbouring vertex collapses and this vertex needs to be updated
	// indexCollapse index of collapsing vertex, indexToCollapseTo index of HEC target, isIndexToCollapseTo if this vertex is the target of the HEC
	void incommingCollapse(int indexCollapse, int indexToCollapseTo, bool isIndexToCollapseTo);

	// Is only called with the vertex that is collapsed to. Updates the triangles that are connected to the collapse.
	void collapseUpdateTriangles(int indexCollapse, int indexToCollapseTo);
private:
	glm::vec3 m_pos = glm::vec3(0);
	bool m_active = false; // true if the vertex is in the mesh false if it should be deleted
	int m_vertexIndex = -1;
	QuadricErrorMetric m_error;
};

inline int Vertex::getVertexIndex() const
{
	return m_vertexIndex;
}
inline glm::vec3 Vertex::getPos() const {
	return m_pos;
}
