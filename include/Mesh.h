/*************************************************
*	Mesh.h
*
*	Release: July 2011
*	Update: July 2022
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
*************************************************/

#pragma once
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <float.h>

//#include <gl/glut.h>

class Vector
{
public:
	Vector(void);
	Vector(const float *v);
	Vector(const double *v);
	Vector(const float *v1, const float *v2);
	Vector(const float v1, const float v2, const float v3);
	Vector(const int v1, const int v2, const int v3);
	~Vector(void);
	Vector cross(const Vector v);
	float norm(void);
	Vector unit(void);
	Vector operator +(const Vector &v);
	Vector operator -(const Vector &v);
	Vector operator -(void);
	Vector operator *(const float v);
	float operator *(const Vector &v);
	Vector operator /(const float v);
	Vector operator =(const float *v);
	Vector operator =(const float v);
	bool operator ==(const Vector &v) const;
	bool operator !=(const Vector &v) const;
	bool operator <(const Vector &v) const;
	bool operator >(const Vector &v) const;
	float operator [](const int id);
	const Vector &operator +=(const Vector &v);
	const Vector &operator -=(const Vector &v);
	const Vector &operator *=(const float v);
	const Vector &operator /=(const float v);
	const float *fv(void);

private:
	float m_vector[3];
};

class Vertex
{
public:
	Vertex(void);
	Vertex(const int id, float *array);
	~Vertex(void);
	void setVertex(const float *v);
	void setList(const int *list, const int n);
	float operator[] (const int id) const;
	const int *list(void) const;
	int list(const int index) const;
	int id(void) const;
	int nNeighbor(void) const;
	const float *fv(void) const;

private:
	float *m_vertex;
	int m_id;
	int *m_list;
	int m_nNeighbor;
	bool m_pre_alloc;
};

class Normal
{
public:
	Normal(void);
	Normal(float *array);
	Normal(const int id, float *array);
	Normal(const float *v);
	~Normal(void);
	void setNormal(const float *v);
	void normalize(void);
	const float *fv(void) const;
	const Normal &operator +=(const Normal &n);
	const Normal &operator *=(const float s);
	float operator[] (const int id) const;

private:
	float *m_normal;
	int m_id;
	bool m_pre_alloc;
};

class Face
{
public:
	Face(void);
	Face(const int id, int *array1, float *array2);
	~Face(void);
	void setVertex(const Vertex **v);
	void setNormal(const Normal **v);
	void setList(const int *list);
	void updateFaceNormal(void);
	int id(void);
	const int *list(void) const;
	int list(const int index) const;
	const Vertex *vertex(const int index) const;
	const Normal *normal(const int index) const;
	const Normal &faceNormal(void) const;
	int operator[] (const int id) const;

private:
	Vertex *m_vertex[3];
	Normal *m_normal[3];
	Normal m_face_normal;
	int *m_list;
	int m_id;
	bool m_pre_alloc;
};

class Mesh
{
public:
	Mesh(void);
	~Mesh(void);
	void draw(void);
	void openFile(const char *filename);
	void saveFile(const char *filename, const char *format, bool normal = false);
	void scaling(const float factor);
	void updateNormal(void);
	void centering(void);
	void unit(void);
	void rotation(const float *axis, float theta);
	void rotation(const float *axis, float theta, float *v);
	void setMesh(const float *vertex, const int *face, const float *normal, int nVertex, int nFace, int nNormal, bool hasNormal);
	int nFace(void) const;
	int nVertex(void) const;
	const Face **face(void) const;
	const Face *face(const int index) const;
	const Vertex **vertex(void) const;
	const Vertex *vertex(const int index) const;
	const Normal **normal(void) const;
	const Normal *normal(const int index) const;

private:
	void connectivity(void);

private:
	int m_nVertex;
	int m_nFace;
	int m_nNormal;
	Vertex **m_vertex;
	Normal **m_normal;
	Face **m_face;
	float *m_vertex_array;
	float *m_normal_array;
	float *m_face_normal_array;
	int *m_face_array;
};
