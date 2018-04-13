/*************************************************
*	MeshIO.cpp
*
*	Release: July 2011
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include "MeshIO.h"

MeshIO::MeshIO(void)
{
	m_nVertex = 0;
	m_nFace = 0;
	m_nNormal = 0;
	m_hasNormal = false;
	m_hasColor = false;
	m_color = NULL;
}

MeshIO::~MeshIO(void)
{
	delete [] m_vertex;
	delete [] m_normal;
	delete [] m_face;
}

int MeshIO::nVertex(void)
{
	return m_nVertex;
}

int MeshIO::nNormal(void)
{
	return m_nNormal;
}

int MeshIO::nFace(void)
{
	return m_nFace;
}

const float *MeshIO::vertex(int id)
{
	return (const float *)&m_vertex[id * 3];
}

const float *MeshIO::normal(int id)
{
	return (const float *)&m_normal[id * 3];
}

const float *MeshIO::color(int id)
{
	return (const float *)&m_color[id * 3];
}

const int *MeshIO::face(int id)
{
	if (m_hasNormal) return (const int *)&m_face[id * 6];
	else return (const int *)&m_face[id * 3];
}

void MeshIO::read(const char *filename)
{
	return;
}

bool MeshIO::hasNormal(void)
{
	return m_hasNormal;
}

bool MeshIO::hasColor(void)
{
	return m_hasColor;
}

void MeshIO::initArray(void)
{
	m_vertex = new float[m_nVertex * 3];
	m_normal = new float[m_nNormal * 3];
	if (m_hasNormal) m_face = new int[m_nFace * 6];
	else m_face = new int[m_nFace * 3];
	if (m_hasColor) m_color = new float[m_nVertex * 3];
}
