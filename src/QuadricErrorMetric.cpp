#include "QuadricErrorMetric.h"
#include "Vertex.h"
#include "ParallelMesh.h"

#undef max

QuadricErrorMetric::QuadricErrorMetric() { 
	m_q = glm::dmat4x4(0);
	m_index = -1;
}

void QuadricErrorMetric::init(Vertex& v)
{
	auto& parent = ParallelMesh::getInstance();

	auto& faceIndices = parent.m_faceIndices;
	auto& adjacentVerticesIndices = parent.m_adjacentVerticesIndices;
#ifdef TRIANGLE_INDICES
	auto& faces = parent.m_faces;
#endif

	{
#ifdef TRIANGLE_INDICES
		int startPos = (slotSize + 1) * v.getVertexIndex();
#else

		int startPos = (3 * slotSize + 1) * v.getVertexIndex();
#endif
		int start = faceIndices[startPos];
		int currentPos = startPos + 1;
		int ctr = 0;
		bool stop = false;

		glm::dvec3 v0, v1, v2;
		glm::dvec3 faceNormal;
		if (faceIndices[currentPos] == -1) { // No Faces Remove This
			m_error = 0;
			m_index = v.getVertexIndex();
		}
		while (!stop && faceIndices[currentPos] != -1) {
#ifndef NO_DEBUG_CHECKS
			if (faceIndices[currentPos] == faceIndices[currentPos + 1] || faceIndices[currentPos] == faceIndices[currentPos + 2] || faceIndices[currentPos + 1] == faceIndices[currentPos + 2]) {
				volatile int i1 = faceIndices[currentPos];
				volatile int i2 = faceIndices[currentPos + 1];
				volatile int i3 = faceIndices[currentPos + 2];
				__debugbreak();
			}
#endif
#ifndef TRIANGLE_INDICES
			v0 = parent.getVertex(faceIndices[currentPos]).getPos();
			v1 = parent.getVertex(faceIndices[currentPos + 1]).getPos();
			v2 = parent.getVertex(faceIndices[currentPos + 2]).getPos();
#else
			v0 = parent.getVertex(faces[faceIndices[currentPos]].x).getPos();
			v1 = parent.getVertex(faces[faceIndices[currentPos]].y).getPos();
			v2 = parent.getVertex(faces[faceIndices[currentPos]].z).getPos();
#endif

			faceNormal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

			if (faceNormal.x == faceNormal.x) //nan happens f.e. for doubled vertices
			{

				double a, b, c, d;
				a = faceNormal.x;
				b = faceNormal.y;
				c = faceNormal.z;

				d = -(a * v0.x + b * v0.y + c * v0.z);


				/*m_11 += a * a;
				m_12 += a * b;
				m_13 += a * c;
				m_14 += a * d;
				m_22 += b * b;
				m_23 += b * c;
				m_24 += b * d;
				m_33 += c * c;
				m_34 += c * d;
				m_44 += d * d;*/

				glm::dmat4x4 k (
					a * a, a * b, a * c, a * d,
					b * a, b * b, b * c, b * d,
					c * a, c * b, c * c, c * d,
					d * a, d * b, d * c, d * d
				);			

				m_q += k;
			}
#ifndef NO_DEBUG_CHECKS
			else {
				//__debugbreak();
			}
#endif
			if (++ctr != slotSize) {
#ifndef TRIANGLE_INDICES
				currentPos += 3;
#else
				currentPos++;
#endif
			}
			else {
				ctr = 0;
				if (start != -1)
				{
					startPos = start;
					start = faceIndices[startPos];
					currentPos = startPos + 1;
				}
				else {
					stop = true;
				}
			}
		}
	}
	m_error = std::numeric_limits<float>::max();
	{
		int startPos = (slotSize + 1) * v.getVertexIndex();
		int start = adjacentVerticesIndices[startPos];
		int currentPos = startPos + 1;
		int ctr = 0;
		bool stop = false;

		//glm::dmat4x4 m_q = getQuadric();
		while (!stop && adjacentVerticesIndices[currentPos] != -1) {

			if (v.getVertexIndex() != adjacentVerticesIndices[currentPos]) // we dont want to calculate our own error (0)
			{
				int vId = adjacentVerticesIndices[currentPos];
				glm::dvec3 pos = parent.getVertexPos(adjacentVerticesIndices[currentPos]);

				glm::dvec4 vq = glm::dvec4(pos, 1) * m_q;
				float error = static_cast<float>(vq.x * pos.x + vq.y * pos.y + vq.z * pos.z + vq.w);
				if (error < 0.f)
				{
					// Ignore negative errors that happen because of precision errors 
					//__debugbreak();
				}
				if (error < m_error)
				{
					m_error = error;
					m_index = adjacentVerticesIndices[currentPos];
				}
			}



			if (++ctr != slotSize) {
				currentPos += 1;
			}
			else {
				ctr = 0;
				if (start != -1)
				{
					startPos = start;
					start = adjacentVerticesIndices[startPos];
					currentPos = startPos + 1;
				}
				else {
					stop = true;
				}
			}
		}
	}
}
void QuadricErrorMetric::updateQuadric(glm::dmat4x4 q2) {
	m_q += q2;
	//glm::dmat4x4 quadricb = v2.getErrorMetric().getQuadric();
	/*m_11 += q2[0][0];
	m_12 += q2[0][1];
	m_13 += q2[0][2];
	m_14 += q2[0][3];
	m_22 += q2[1][1];
	m_23 += q2[1][2];
	m_24 += q2[1][3];
	m_33 += q2[2][2];
	m_34 += q2[2][3];
	m_44 += q2[3][3];*/
	/*m_q[0] += q2[0][0];
	m_q[1] += q2[0][1];
	m_q[2] += q2[0][2];
	m_q[3] += q2[0][3];
	m_q[4] += q2[1][1];
	m_q[5] += q2[1][2];
	m_q[6] += q2[1][3];
	m_q[7] += q2[2][2];
	m_q[8] += q2[2][3];
	m_q[9] += q2[3][3];*/
}
float QuadricErrorMetric::setNextBestError(Vertex& v)

{

	bool once = false;
	auto& parent = ParallelMesh::getInstance();
	auto& adjacentVerticesIndices = parent.m_adjacentVerticesIndices;
	auto& faceIndices = parent.m_faceIndices;
#ifdef TRIANGLE_INDICES
	auto& faces = parent.m_faces;
#endif // !TRIANGLE_INDICES
	{
		float oldError = m_error;
		int oldIndex = m_index;
		m_error = std::numeric_limits<float>::max();
		{
#ifndef TRIANGLE_INDICES
			int startPos = (3 * slotSize + 1) * v.getVertexIndex();
#else
			int startPos = (slotSize + 1) * v.getVertexIndex();
#endif
			int currentPos = startPos + 1;
			if (faceIndices[currentPos] == -1) { // No Faces Remove This
				m_error = 0;
				m_index = v.getVertexIndex();
			}
		}
		{
			int startPos = (slotSize + 1) * v.getVertexIndex();
			int start = adjacentVerticesIndices[startPos];
			int currentPos = startPos + 1;
			int ctr = 0;
			bool stop = false;
			//glm::dmat4x4 m_q = getQuadric();

			while (!stop && adjacentVerticesIndices[currentPos] != -1) {
				if (v.getVertexIndex() != adjacentVerticesIndices[currentPos]) // we dont want to calculate our own error (0)
				{
					glm::dvec3 pos = parent.getVertex(adjacentVerticesIndices[currentPos]).getPos();
					glm::dvec4 vq = glm::dvec4(pos, 1) * (m_q + parent.getVertex(adjacentVerticesIndices[currentPos]).getErrorMetric().getQuadric());
					float error = static_cast<float>(vq.x * pos.x + vq.y * pos.y + vq.z * pos.z + vq.w);
					if (error < 0.f)
					{
						// Ignore negative errors that happen because of precision errors 
						//__debugbreak();
					}
					if (error < m_error)
					{
						if (error > oldError)
						{
							m_error = error;
							m_index = adjacentVerticesIndices[currentPos];
							once = true;
						}
					}
				}


				if (++ctr != slotSize) {
					currentPos += 1;
				}
				else {
					ctr = 0;
					if (start != -1)
					{
						startPos = start;
						start = adjacentVerticesIndices[startPos];
						currentPos = startPos + 1;
					}
					else {
						stop = true;
					}
				}
			}
		}
	}
	if (once)
		return m_error;
	else
		return 1000000000000000000000000.0f;
}

void QuadricErrorMetric::update(Vertex& v)
{
		

	auto& parent = ParallelMesh::getInstance();
	auto& adjacentVerticesIndices = parent.m_adjacentVerticesIndices;
	auto& faceIndices = parent.m_faceIndices;
#ifdef TRIANGLE_INDICES
	auto& faces = parent.m_faces;
#endif // #ifdef TRIANGLE_INDICES
	{
		m_error = std::numeric_limits<float>::max();
		{
#ifndef TRIANGLE_INDICES
			int startPos = (3 * slotSize + 1) * v.getVertexIndex();
#else
			int startPos = (slotSize + 1) * v.getVertexIndex();
#endif
			int currentPos = startPos + 1;
			if (faceIndices[currentPos] == -1) { // No Faces Remove This
				m_error = 0;
				m_index = v.getVertexIndex();
			}
		}
		{
			int startPos = (slotSize + 1) * v.getVertexIndex();
			int start = adjacentVerticesIndices[startPos];
			int currentPos = startPos + 1;
			int ctr = 0;
			bool stop = false;
			//glm::dmat4x4 m_q = getQuadric();

			while (!stop && adjacentVerticesIndices[currentPos] != -1) {
				if (v.getVertexIndex() != adjacentVerticesIndices[currentPos]) // we dont want to calculate our own error (0)
				{
					glm::dvec3 pos = parent.getVertex(adjacentVerticesIndices[currentPos]).getPos();
					glm::dvec4 vq = glm::dvec4(pos, 1) * (m_q + parent.getVertex(adjacentVerticesIndices[currentPos]).getErrorMetric().getQuadric());
					float error = static_cast<float>(vq.x * pos.x + vq.y * pos.y + vq.z * pos.z + vq.w);
					if (error < 0.f)
					{
						// Ignore negative errors that happen because of precision errors 
						//__debugbreak();
					}
					if (error < m_error)
					{
						m_error = error;
						m_index = adjacentVerticesIndices[currentPos];
					}
				}


				if (++ctr != slotSize) {
					currentPos += 1;
				}
				else {
					ctr = 0;
					if (start != -1)
					{
						startPos = start;
						start = adjacentVerticesIndices[startPos];
						currentPos = startPos + 1;
					}
					else {
						stop = true;
					}
				}
			}
		}
	}

#ifndef NO_DEBUG_CHECKS
	if (m_index == -1) {
		{
			__debugbreak();
			volatile int startPos = (slotSize + 1) * v.getVertexIndex();
			volatile int start = adjacentVerticesIndices[startPos];
			volatile int currentPos = startPos + 1;
			volatile int ctr = 0;
			volatile bool stop = false;
			volatile int qq = v.getVertexIndex();
			//glm::dmat4x4 m_q = getQuadric();
			while (!stop && adjacentVerticesIndices[currentPos] != -1) {
				volatile float test = adjacentVerticesIndices[currentPos];
				if (v.getVertexIndex() != adjacentVerticesIndices[currentPos]) // we dont want to calculate our own error (0)
				{
					glm::dvec3 pos = parent.getVertex(adjacentVerticesIndices[currentPos]).getPos();
					glm::dvec4 vq = glm::dvec4(pos, 1) * (m_q + parent.getVertex(adjacentVerticesIndices[currentPos]).getErrorMetric().getQuadric());
					volatile float x1 = vq.x;
					volatile float x2 = vq.y;
					volatile float x3 = vq.z;
					volatile float error = vq.x * pos.x + vq.y * pos.y + vq.z * pos.z + vq.w;
					if (error < 0.f)
					{
						// Ignore negative errors that happen because of precision errors 
						//__debugbreak();
					}
					if (error < m_error)
					{
						m_error = error;
						m_index = adjacentVerticesIndices[currentPos];
					}
				}



				if (++ctr != slotSize) {
					currentPos += 1;
				}
				else {
					ctr = 0;
					if (start != -1)
					{
						startPos = start;
						start = adjacentVerticesIndices[startPos];
						currentPos = startPos + 1;
					}
					else {
						stop = true;
					}
				}
			}
		}
	}
#endif
}

void QuadricErrorMetric::updateSingle(int newVertex)
{


	auto& parent = ParallelMesh::getInstance();

	{
		glm::dvec3 pos = parent.getVertex(newVertex).getPos();

		glm::dvec4 vq = glm::dvec4(pos, 1) * (m_q + parent.getVertex(newVertex).getErrorMetric().getQuadric());
		//glm::dvec4 vq = glm::dvec4(pos, 1) * getQuadric();
		float error = static_cast<float>(vq.x * pos.x + vq.y * pos.y + vq.z * pos.z + vq.w);
		if (error < 0.f)
		{
			// Ignore negative errors that happen because of precision errors 
			//__debugbreak();
		}
		if (error < m_error)
		{
			m_error = error;
			m_index = newVertex;
		}
	}
}