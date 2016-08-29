/*************************************************
*	MeshIO.h
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
#include <iostream>
#include <fstream>
#include <string>
#include "Mesh.h"

class MeshIO
{
public:
	MeshIO(void);
	~MeshIO(void);
	virtual void read(const char *filename);
	void save(const char *filename, const Mesh *mesh);
	int nVertex(void);
	int nNormal(void);
	int nFace(void);
	bool hasNormal(void);
	bool hasColor(void);
	const float *vertex(int id);
	const float *normal(int id);
	const float *color(int id);
	const int *face(int id);

protected:
	void initArray(void);

protected:
	int m_nVertex;
	int m_nFace;
	int m_nNormal;
	float *m_vertex;
	float *m_color;
	float *m_normal;
	int *m_face;
	bool m_hasNormal;
	bool m_hasColor;
};

