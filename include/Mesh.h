/*************************************************
*	Mesh.h
*
*	Release: July 2011
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#pragma once
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
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
	Vertex(const int id);
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
	float m_vertex[3];
	int m_id;
	int *m_list;
	int m_nNeighbor;
};

class Normal
{
public:
	Normal(void);
	Normal(const int id);
	Normal(const float *v);
	~Normal(void);
	void setNormal(const float *v);
	void normalize(void);
	const float *fv(void) const;
	const Normal &operator +=(const Normal &n);
	const Normal &operator *=(const float s);
	float operator[] (const int id) const;

private:
	float m_normal[3];
	int m_id;
};

class Face
{
public:
	Face(void);
	Face(const int id);
	~Face(void);
	void setVertex(const Vertex **v);
	void setNormal(const Normal **v);
	void setList(const int *list);
	int id(void);
	const int *list(void) const;
	int list(const int index) const;
	const Vertex *vertex(const int index) const;
	const Normal *normal(const int index) const;
	const Normal faceNormal(void) const;
	int operator[] (const int id) const;

private:
	Vertex *m_vertex[3];
	Normal *m_normal[3];
	int m_list[3];
	int m_id;
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
	void rotation(const float *axis, float theta);
	void rotation(const float *axis, float theta, float *v);
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
};
