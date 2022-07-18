/*************************************************
*	Mesh.cpp
*
*	Release: July 2011
*	Update: July 2022
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
*************************************************/

#include "Mesh.h"
#include "MNIObjIO.h"
#include "VtkIO.h"
#include "string.h"

using namespace std;

Vector::Vector(void)
{
	m_vector[0] = m_vector[1] = m_vector[2] = 0;
}

Vector::Vector(const float *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
}

Vector::Vector(const double *v)
{
	m_vector[0] = (float)v[0];
	m_vector[1] = (float)v[1];
	m_vector[2] = (float)v[2];
}

Vector::Vector(const float *v1, const float *v2)
{
	m_vector[0] = v2[0] - v1[0];
	m_vector[1] = v2[1] - v1[1];
	m_vector[2] = v2[2] - v1[2];
}

Vector::Vector(const float v1, const float v2, const float v3)
{
	m_vector[0] = v1;
	m_vector[1] = v2;
	m_vector[2] = v3;
}

Vector::Vector(const int v1, const int v2, const int v3)
{
	m_vector[0] = (const float)v1;
	m_vector[1] = (const float)v2;
	m_vector[2] = (const float)v3;
}

Vector::~Vector(void)
{
}

Vector Vector::cross(const Vector v)
{
	Vector p;
	p.m_vector[0] = m_vector[1] * v.m_vector[2] - m_vector[2] * v.m_vector[1];
	p.m_vector[1] = m_vector[2] * v.m_vector[0] - m_vector[0] * v.m_vector[2];
	p.m_vector[2] = m_vector[0] * v.m_vector[1] - m_vector[1] * v.m_vector[0];
	return p;
}

float Vector::norm(void)
{
	return sqrt(m_vector[0] * m_vector[0] + m_vector[1] * m_vector[1] + m_vector[2] * m_vector[2]);
}

Vector Vector::unit(void)
{
	float s = norm();
	if (s > 0)
	{
		m_vector[0] /= s;
		m_vector[1] /= s;
		m_vector[2] /= s;
	}
	return *this;
}

const float * Vector::fv(void)
{
	return m_vector;
}

Vector Vector::operator +(const Vector &v)
{
	return Vector(m_vector[0] + v.m_vector[0], m_vector[1] + v.m_vector[1], m_vector[2] + v.m_vector[2]);
}

Vector Vector::operator -(const Vector &v)
{
	return Vector(m_vector[0] - v.m_vector[0], m_vector[1] - v.m_vector[1], m_vector[2] - v.m_vector[2]);
}

Vector Vector::operator -(void)
{
	return Vector(-m_vector[0], -m_vector[1], -m_vector[2]);
}

Vector Vector::operator *(const float v)
{
	return Vector(m_vector[0] * v, m_vector[1] * v, m_vector[2] * v);
}

float Vector::operator *(const Vector &v)
{
	return m_vector[0] * v.m_vector[0] + m_vector[1] * v.m_vector[1] + m_vector[2] * v.m_vector[2];
}

Vector Vector::operator /(const float v)
{
	return Vector(m_vector[0] / v, m_vector[1] / v, m_vector[2] / v);
}

Vector Vector::operator =(const float *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
	return *this;
}

Vector Vector::operator =(const float v)
{
	m_vector[0] = v;
	m_vector[1] = 0;
	m_vector[2] = 0;
	return *this;
}

const Vector & Vector::operator +=(const Vector &v)
{
	m_vector[0] += v.m_vector[0];
	m_vector[1] += v.m_vector[1];
	m_vector[2] += v.m_vector[2];
	return *this;
}

const Vector & Vector::operator -=(const Vector &v)
{
	m_vector[0] -= v.m_vector[0];
	m_vector[1] -= v.m_vector[1];
	m_vector[2] -= v.m_vector[2];
	return *this;
}

const Vector & Vector::operator *=(const float v)
{
	m_vector[0] *= v; m_vector[1] *= v; m_vector[2] *= v;
	return *this;
}

const Vector & Vector::operator /=(const float v)
{
	m_vector[0] /= v; m_vector[1] /= v; m_vector[2] /= v;
	return *this;
}

bool Vector::operator ==(const Vector &v) const
{
	return (m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

float Vector::operator[] (const int id)
{
	return m_vector[id];
}

bool Vector::operator !=(const Vector &v) const
{
	return !(m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

bool Vector::operator <(const Vector &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] < v.m_vector[2];
		else
			return m_vector[1] < v.m_vector[1];
	else
		return m_vector[0] < v.m_vector[0];
}

bool Vector::operator >(const Vector &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] > v.m_vector[2];
		else
			return m_vector[1] > v.m_vector[1];
	else
		return m_vector[0] > v.m_vector[0];
}

Vertex::Vertex(void)
{
	m_pre_alloc = false;
	m_vertex = new float[3];
	m_vertex[0] = m_vertex[1] = m_vertex[2] = 0;
	m_id = 0;
	m_nNeighbor = 0;
	m_list = NULL;
}

Vertex::Vertex(const int id, float *array)
{
	m_pre_alloc = true;
	m_id = id;
	m_nNeighbor = 0;
	m_list = NULL;
	m_vertex = array;
}

Vertex::~Vertex(void)
{
	delete [] m_list;
	if (!m_pre_alloc) delete [] m_vertex;
}

const float * Vertex::fv(void) const
{
	return m_vertex;
}

float Vertex::operator[] (const int id) const
{
	return m_vertex[id];
}

void Vertex::setVertex(const float *v)
{
	m_vertex[0] = v[0]; m_vertex[1] = v[1]; m_vertex[2] = v[2];
}

void Vertex::setList(const int *list, const int n)
{
	m_nNeighbor = n;
	if (m_list != NULL) delete [] m_list;
	m_list = new int[n];
	memcpy(m_list, list, sizeof(int) * n);
}

void Vertex::setPreAlloc(float *array)
{
	if (!m_pre_alloc) delete [] m_vertex;
	m_pre_alloc = true;
	m_vertex = array;
}

int Vertex::id(void) const
{
	return m_id;
}

const int * Vertex::list(void) const
{
	return m_list;
}

int Vertex::list(const int index) const
{
	return m_list[index];
}

int Vertex::nNeighbor(void) const
{
	return m_nNeighbor;
}

Normal::Normal(void)
{
	m_pre_alloc = false;
	m_normal = new float[3];
	m_normal[0] = m_normal[1] = m_normal[2] = 0;
	m_id = 0;
}

Normal::Normal(float *array)
{
	m_pre_alloc = true;
	m_normal = array;
	m_id = 0;
}

Normal::Normal(const int id, float *array)
{
	m_pre_alloc = true;
	m_normal = array;
	m_normal[0] = m_normal[1] = m_normal[2] = 0;
	m_id = id;
}

Normal::Normal(const float *v)
{
	m_pre_alloc = false;
	m_normal = new float[3];
	setNormal(v);
}

Normal::~Normal(void)
{
	if (!m_pre_alloc) delete [] m_normal;
}

const float * Normal::fv(void) const
{
	return m_normal;
}

void Normal::normalize(void)
{
	Vector v(m_normal);
	v.unit();
	setNormal(v.fv());
}

void Normal::setNormal(const float *v)
{
	m_normal[0] = v[0]; m_normal[1] = v[1]; m_normal[2] = v[2];
}

void Normal::setPreAlloc(float *array)
{
	if (!m_pre_alloc) delete [] m_normal;
	m_pre_alloc = true;
	m_normal = array;
}

const Normal & Normal::operator +=(const Normal &n)
{
	m_normal[0] += n.m_normal[0];
	m_normal[1] += n.m_normal[1];
	m_normal[2] += n.m_normal[2];

	return *this;
}

const Normal & Normal::operator *=(const float s)
{
	m_normal[0] *= s;
	m_normal[1] *= s;
	m_normal[2] *= s;

	return *this;
}

float Normal::operator[] (const int id) const
{
	return m_normal[id];
}

Face::Face(void)
{
	m_pre_alloc = false;
	m_list = new int[3];
	m_id = 0;
}

Face::Face(const int id, int *array1, float *array2)
{
	m_pre_alloc = true;
	m_list = array1;
	m_face_normal.setPreAlloc(array2);
	m_id = id;
}

Face::~Face(void)
{
	if (!m_pre_alloc) delete [] m_list;
}

void Face::setVertex(const Vertex **v)
{
	m_vertex[0] = (Vertex *)v[0]; m_vertex[1] = (Vertex *)v[1]; m_vertex[2] = (Vertex *)v[2];
	updateFaceNormal();
}

void Face::setPreAlloc(int *array1, float *array2)
{
	if (!m_pre_alloc) delete [] m_list;
	m_pre_alloc = true;
	m_list = array1;
	m_face_normal.setPreAlloc(array2);
}

const Vertex * Face::vertex(const int index) const
{
	return m_vertex[index];
}

void Face::setNormal(const Normal **v)
{
	m_normal[0] = (Normal *)v[0]; m_normal[1] = (Normal *)v[1]; m_normal[2] = (Normal *)v[2];
}

void Face::setList(const int *list)
{
	m_list[0] = list[0];
	m_list[1] = list[1];
	m_list[2] = list[2];
}

void Face::updateFaceNormal(void)
{
	Vector v1(m_vertex[1]->fv(), m_vertex[0]->fv());
	Vector v2(m_vertex[2]->fv(), m_vertex[1]->fv());
	Vector v3 = v1.cross(v2);
	v3.unit();
	m_face_normal.setNormal(v3.fv());
}
const int *Face::list(void) const
{
	return (const int *)m_list;
}

int Face::list(const int index) const
{
	return (const int)m_list[index];
}

int Face::operator[] (const int id) const
{
	return m_list[id];
}

const Normal * Face::normal(const int index) const
{
	return m_normal[index];
}

const Normal & Face::faceNormal(void) const
{
	return m_face_normal;
}

int Face::id(void)
{
	return m_id;
}

Mesh::Mesh(void)
{
	m_nVertex = 0;
	m_nFace = 0;
	m_nNormal = 0;

	m_vertex = NULL;
	m_normal = NULL;
	m_face = NULL;

	m_vertex_array = NULL;
	m_normal_array = NULL;
	m_face_normal_array = NULL;
	m_face_array = NULL;
}

Mesh::~Mesh(void)
{
	if (m_vertex != NULL)
	{
		for (int i = 0; i < m_nVertex; i++) delete m_vertex[i];
		delete [] m_vertex;
	}
	if (m_normal != NULL)
	{
		for (int i = 0; i < m_nNormal; i++) delete m_normal[i];
		delete [] m_normal;
	}
	if (m_face != NULL)
	{
		for (int i = 0; i < m_nFace; i++) delete m_face[i];
		delete [] m_face;
	}
	if (m_vertex_array != NULL) delete [] m_vertex_array;
	if (m_normal_array != NULL) delete [] m_normal_array;
	if (m_face_normal_array != NULL) delete [] m_face_normal_array;
	if (m_face_array != NULL) delete [] m_face_array;
}

const Vertex ** Mesh::vertex(void) const
{
	return (const Vertex **)m_vertex;
}

const Vertex * Mesh::vertex(int index) const
{
	return (const Vertex *)m_vertex[index];
}

const Normal ** Mesh::normal(void) const
{
	return (const Normal **)m_normal;
}

const Normal * Mesh::normal(int index) const
{
	return (const Normal *)m_normal[index];
}

const Face ** Mesh::face(void) const
{
	return (const Face **)m_face;
}

const Face * Mesh::face(int index) const
{
	return (const Face *)m_face[index];
}

int Mesh::nFace(void) const
{
	return m_nFace;
}

int Mesh::nVertex(void) const
{
	return m_nVertex;
}

void Mesh::rotation(const float *axis, float theta)
{
	if (theta == 0) return;
	float P[9] = {axis[0] * axis[0], axis[0] * axis[1], axis[0] * axis[2],
			axis[0] * axis[1], axis[1] * axis[1], axis[1] * axis[2],
			axis[0] * axis[2], axis[1] * axis[2], axis[2] * axis[2]};
	float Q[9] = {0, -axis[2], axis[1], axis[2], 0, -axis[0], -axis[1], axis[0], 0};
	float I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	float c = cos(theta);
	float s = sin(theta);
	float mat[9];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			mat[i * 3 + j] = P[i * 3 + j] + (I[i * 3 + j] - P[i * 3 + j]) * c + Q[i * 3 + j] * s;

	for (int n = 0; n < m_nVertex; n++)
	{
		float v[3] = {0, 0, 0};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				v[i] += mat[i * 3 + j] * m_vertex[n]->fv()[j];
		m_vertex[n]->setVertex(v);
	}
}

void Mesh::centering(void)
{
	float center[3] = {0, 0, 0};
	for (int n = 0; n < m_nVertex; n++)
	{
		center[0] += m_vertex[n]->fv()[0];
		center[1] += m_vertex[n]->fv()[1];
		center[2] += m_vertex[n]->fv()[2];
	}
	center[0] /= m_nVertex; center[1] /= m_nVertex; center[2] /= m_nVertex;
	for (int n = 0; n < m_nVertex; n++)
	{
		float v[3] = {m_vertex[n]->fv()[0] - center[0], m_vertex[n]->fv()[1] - center[1], m_vertex[n]->fv()[2] - center[2]};
		m_vertex[n]->setVertex(v);
	}
}

void Mesh::scaling(const float factor)
{
	if (factor == 0) return;
	for (int n = 0; n < m_nVertex; n++)
	{
		float v[3] = {m_vertex[n]->fv()[0] * factor, m_vertex[n]->fv()[1] * factor, m_vertex[n]->fv()[2] * factor};
		m_vertex[n]->setVertex(v);
	}
}

void Mesh::rotation(const float *axis, float theta, float *v)
{
	if (theta == 0) return;
	float P[9] = {axis[0] * axis[0], axis[0] * axis[1], axis[0] * axis[2],
			axis[0] * axis[1], axis[1] * axis[1], axis[1] * axis[2],
			axis[0] * axis[2], axis[1] * axis[2], axis[2] * axis[2]};
	float Q[9] = {0, -axis[2], axis[1], axis[2], 0, -axis[0], -axis[1], axis[0], 0};
	float I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	float c = cos(theta);
	float s = sin(theta);
	float fv[3] = {v[0], v[1], v[2]};
	v[0] = v[1] = v[2] = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			v[i] += (P[i * 3 + j] + (I[i * 3 + j] - P[i * 3 + j]) * c + Q[i * 3 + j] * s) * fv[j];
}

void Mesh::setMesh(const float *vertex, const int *face, const float *normal, int nVertex, int nFace, int nNormal, bool hasNormal)
{
	if (m_vertex != NULL)
	{
		for (int i = 0; i < m_nVertex; i++) delete m_vertex[i];
		delete [] m_vertex;
	}
	if (m_normal != NULL)
	{
		for (int i = 0; i < m_nNormal; i++) delete m_normal[i];
		delete [] m_normal;
	}
	if (m_face != NULL)
	{
		for (int i = 0; i < m_nFace; i++) delete m_face[i];
		delete [] m_face;
	}
	if (m_vertex_array != NULL) delete [] m_vertex_array;
	if (m_normal_array != NULL) delete [] m_normal_array;
	if (m_face_normal_array != NULL) delete [] m_face_normal_array;
	if (m_face_array != NULL) delete [] m_face_array;

	m_nVertex = nVertex;
	m_nFace = nFace;
	m_nNormal = (hasNormal) ? nNormal: nVertex;

	//cout << "# of vertices: " << m_nVertex << endl;
	//cout << "# of faces: " << m_nFace << endl;

	if (m_nNormal != m_nVertex) m_nNormal = m_nVertex;
	m_vertex = new Vertex*[m_nVertex];
	m_normal = new Normal*[m_nNormal];
	m_face = new Face*[m_nFace];

	m_vertex_array = new float[m_nVertex * 3];
	m_normal_array = new float[m_nNormal * 3];
	m_face_normal_array = new float[m_nFace * 3];
	m_face_array = new int[m_nFace * 3];

	for (int i = 0; i < m_nVertex; i++) m_vertex[i] = new Vertex(i, &m_vertex_array[i * 3]);
	for (int i = 0; i < m_nNormal; i++) m_normal[i] = new Normal(i, &m_normal_array[i * 3]);
	for (int i = 0; i < m_nFace; i++) m_face[i] = new Face(i, &m_face_array[i * 3], &m_face_normal_array[i * 3]);

	for (int i = 0; i < m_nVertex; i++)
	{
		const float *v = &vertex[i * 3];
		m_vertex[i]->setVertex(v);
	}
	if (hasNormal)
	{
		for (int i = 0; i < m_nNormal; i++)
		{
			const float *v = &normal[i * 3];
			m_normal[i]->setNormal(v);
		}
	}
	for (int i = 0; i < m_nFace; i++)
	{
		const int *idx = (hasNormal) ? &face[i * 6]: &face[i * 3];
		const Vertex *v[3];
		const Normal *n[3];
		if (hasNormal)
		{
			v[0] = m_vertex[idx[0]]; v[1] = m_vertex[idx[1]]; v[2] = m_vertex[idx[2]];
			n[0] = m_normal[idx[3]]; n[1] = m_normal[idx[4]]; n[2] = m_normal[idx[5]];
			m_face[i]->setVertex(v);
			m_face[i]->setNormal(n);
			m_face[i]->setList(idx);
		}
		else
		{
			v[0] = m_vertex[idx[0]]; v[1] = m_vertex[idx[1]]; v[2] = m_vertex[idx[2]];
			n[0] = m_normal[idx[0]]; n[1] = m_normal[idx[1]]; n[2] = m_normal[idx[2]];
			m_face[i]->setVertex(v);
			m_face[i]->setNormal(n);
			m_face[i]->setList(idx);

			const int *idx = m_face[i]->list();
			Vector e1(m_vertex[idx[1]]->fv(), m_vertex[idx[0]]->fv());
			Vector e2(m_vertex[idx[2]]->fv(), m_vertex[idx[1]]->fv());
			Vector e3(m_vertex[idx[0]]->fv(), m_vertex[idx[2]]->fv());
			Vector FN = e1.cross(e2);
			e1.unit(); e2.unit(); e3.unit();
			float e13 = e1 * -e3;
			if (e13 > 1) e13 = 1;
			if (e13 < -1) e13 = -1;
			float e21 = e2 * -e1;
			if (e21 > 1) e21 = 1;
			if (e21 < -1) e21 = -1;
			float e32 = e3 * -e2;
			if (e32 > 1) e32 = 1;
			if (e32 < -1) e32 = -1;
			float w1 = acos(e13);
			float w2 = acos(e21);
			float w3 = acos(e32);

			*m_normal[idx[0]] += Normal((FN * w1).fv());
			*m_normal[idx[1]] += Normal((FN * w2).fv());
			*m_normal[idx[2]] += Normal((FN * w3).fv());
			/*Normal fn = m_face[i]->faceNormal();
			*m_normal[idx[0]] += fn;
			*m_normal[idx[1]] += fn;
			*m_normal[idx[2]] += fn;*/
		}
	}

	for (int i = 0; i < m_nNormal; i++) m_normal[i]->normalize();

	connectivity();
}

void Mesh::openFile(const char *filename)
{
	ifstream fin(filename);
	//cout << "Filename: " << filename << endl;
	if(fin.fail())
	{
		cout << "Failed to open " << filename << endl;
		return;
	}
	fin.close();

	const char *format = &filename[strlen(filename) - 3];
	MeshIO *mesh;
	if (!strcmp(format, "obj")) mesh = new MNIObjIO(filename);
	else if (!strcmp(format, "vtk")) mesh = new VtkIO(filename);
	else
	{
		cout << "Unsupported format!" << endl;
		return;
	}

	setMesh(mesh->vertex(0), mesh->face(0), mesh->normal(0), mesh->nVertex(), mesh->nFace(), mesh->nNormal(), mesh->hasNormal());
	updateNormal();

	delete mesh;
}

void Mesh::updateNormal(void)
{
	for (int i = 0; i < m_nNormal; i++) *m_normal[i] *= 0;
	for (int i = 0; i < m_nFace; i++)
	{
		const int *idx = m_face[i]->list();
		Vector e1(m_vertex[idx[1]]->fv(), m_vertex[idx[0]]->fv());
		Vector e2(m_vertex[idx[2]]->fv(), m_vertex[idx[1]]->fv());
		Vector e3(m_vertex[idx[0]]->fv(), m_vertex[idx[2]]->fv());
		Vector FN = e1.cross(e2);
		e1.unit(); e2.unit(); e3.unit();
		float e13 = e1 * -e3;
		if (e13 > 1) e13 = 1;
		if (e13 < -1) e13 = -1;
		float e21 = e2 * -e1;
		if (e21 > 1) e21 = 1;
		if (e21 < -1) e21 = -1;
		float e32 = e3 * -e2;
		if (e32 > 1) e32 = 1;
		if (e32 < -1) e32 = -1;
		float w1 = acos(e13);
		float w2 = acos(e21);
		float w3 = acos(e32);

		*m_normal[idx[0]] += Normal((FN * w1).fv());
		*m_normal[idx[1]] += Normal((FN * w2).fv());
		*m_normal[idx[2]] += Normal((FN * w3).fv());
		/*Normal fn = m_face[i]->faceNormal();
		*m_normal[idx[0]] += fn;
		*m_normal[idx[1]] += fn;
		*m_normal[idx[2]] += fn;*/

		m_face[i]->updateFaceNormal();
	}
	for (int i = 0; i < m_nNormal; i++) m_normal[i]->normalize();
}

void Mesh::saveFile(const char *filename, const char *format, bool normal)
{

	if (!strcmp(format, "vtk"))
	{
		VtkIO::save(filename, this, normal);
	}
	/*else if (!strcmp(format, "obj"))
	{
		MNIObjIO::save(filename, this, normal);
	}*/
	else
	{
		cout << "Unsupported format!" << endl;
		return;
	}
}

void Mesh::unit(void)
{
	centering();

	float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
	float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
	for (int n = 0; n < m_nVertex; n++)
	{
		if (minx > m_vertex[n]->fv()[0]) minx = m_vertex[n]->fv()[0];
		if (maxx < m_vertex[n]->fv()[0]) maxx = m_vertex[n]->fv()[0];
		if (miny > m_vertex[n]->fv()[1]) miny = m_vertex[n]->fv()[1];
		if (maxy < m_vertex[n]->fv()[1]) maxy = m_vertex[n]->fv()[1];
		if (minz > m_vertex[n]->fv()[2]) minz = m_vertex[n]->fv()[2];
		if (maxz < m_vertex[n]->fv()[2]) maxz = m_vertex[n]->fv()[2];
	}
	float lenx = 1 / (maxx - minx);
	float leny = 1 / (maxy - miny);
	float lenz = 1 / (maxz - minz);

	float factor;
	factor = (lenx > leny) ? lenx: leny;
	factor = (factor > lenz) ? factor: lenz;

	scaling(factor);
}

void Mesh::connectivity(void)
{
	int *nNeighbor = new int[m_nVertex];
	memset(nNeighbor, 0, sizeof(int) * m_nVertex);

	for (int i = 0; i < m_nFace; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int idx = m_face[i]->vertex(j)->id();
			nNeighbor[idx]++;
		}
	}

	struct link
	{
		int idx;
		int next;
	};
	link **neighbor = new link*[m_nVertex];
	for (int i = 0; i < m_nVertex; i++)
	{
		//cout << i << " " << nNeighbor[i] << endl;
		neighbor[i] = new link[nNeighbor[i]];
	}

	memset(nNeighbor, 0, sizeof(int) * m_nVertex);
	for (int i = 0; i < m_nFace; i++)
	{
		int idx0 = m_face[i]->vertex(0)->id();
		int idx1 = m_face[i]->vertex(1)->id();
		int idx2 = m_face[i]->vertex(2)->id();

		neighbor[idx0][nNeighbor[idx0]].idx = idx1;
		neighbor[idx0][nNeighbor[idx0]].next = idx2;
		nNeighbor[idx0]++;

		neighbor[idx1][nNeighbor[idx1]].idx = idx2;
		neighbor[idx1][nNeighbor[idx1]].next = idx0;
		nNeighbor[idx1]++;

		neighbor[idx2][nNeighbor[idx2]].idx = idx0;
		neighbor[idx2][nNeighbor[idx2]].next = idx1;
		nNeighbor[idx2]++;
	}

	for (int i = 0; i < m_nVertex; i++)
	{
		int *list = new int[nNeighbor[i]];
		int next = 0;
		int j = 0;
		int n = nNeighbor[i];	// pointer list
		bool *visited = new bool[n];
		memset(visited, 0, sizeof(bool) * n);
		while (next < n)
		{
			list[j++] = neighbor[i][next].idx;
			visited[next] = true;

			bool found = false;
			for (int k = 0; k < n; k++)
			{
				if (neighbor[i][k].idx == neighbor[i][next].next)
				{
					next = k;
					found = true;
					break;
				}
			}
			if (!found)
			{
				int *newList = new int[nNeighbor[i] + 1];
				memcpy(newList, list, sizeof(int) * nNeighbor[i]);
				nNeighbor[i]++;
				delete [] list;
				list = newList;
				list[j++] = neighbor[i][next].next;
			}
			if (visited[next])
			{
				for (next = 0; next < n; next++) if (!visited[next]) break;
				if (next == n) break;
			}
		}
		m_vertex[i]->setList(list, nNeighbor[i]);
		delete [] list;
		delete [] visited;
	}

	// free
	delete [] nNeighbor;
	for (int i = 0; i < m_nVertex; i++) delete [] neighbor[i];
	delete [] neighbor;
}
