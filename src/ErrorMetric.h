#pragma once

#include "defines.h"
#include <glm/glm.hpp>

class Vertex;

class ErrorMetric {
public:
	virtual void update(Vertex& v) = 0;
	virtual void updateQuadric(glm::dmat4x4 q2) = 0;
	virtual glm::dmat4x4& getQuadric() = 0;
	virtual void init(Vertex& v) = 0;

	inline float getSmallestError() {
		return m_error;
	}
	inline int getSmallestErrorIndex() {
		return m_index;
	}
protected:
	float m_error;
	int m_index;
};