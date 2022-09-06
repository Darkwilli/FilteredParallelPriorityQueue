#pragma once

#include <glm/glm.hpp>

#include "ErrorMetric.h"

// The QuadricErrorMetric contains the error quadric, the best vertex index to collapse to and its HEC error
class QuadricErrorMetric{
public:
	QuadricErrorMetric();

	void update(Vertex& v);
	void updateSingle(int newVertex);
	void updateQuadric(glm::dmat4x4 q2);
	float setNextBestError(Vertex& v);
	glm::dmat4x4& getQuadric() {
		return m_q;
		/*glm::dmat4x4 quad = {
			{m_11, m_12, m_13, m_14},
			{m_12, m_22, m_23, m_24},
			{m_13, m_23, m_33, m_34},
			{m_14, m_24, m_34, m_44}
		};*/
		/*glm::dmat4 quad (
			m_q[0], m_q[1], m_q[2], m_q[3],
			m_q[1], m_q[4], m_q[5], m_q[6],
			m_q[2], m_q[5], m_q[7], m_q[8],
			m_q[3], m_q[6], m_q[8], m_q[9]
		);*/
		///return quad;
	}
	void init(Vertex& v);
	inline float getSmallestError() {
		return m_error;
	}
	inline int getSmallestErrorIndex() {
		return m_index;
	}
private:
	float m_error = 0;
	int m_index = -1;
	//std::array<double, 10> m_q = { 0,0,0,0,0,0,0,0,0,0 };
	//double m_11, m_12, m_13, m_14, m_22, m_23, m_24, m_33, m_34, m_44;
	//std::unique_ptr<glm::dmat4x4> m_q = std::make_unique<glm::dmat4x4>();
	glm::dmat4x4 m_q;
};