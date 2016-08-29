/*************************************************
*	MNIObjIO.cpp
*
*	Release: March 2015
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include "MNIObjIO.h"

#include "stdlib.h"
#include "iostream"
#include "fstream"
#include "string.h"

using namespace std;

MNIObjIO::MNIObjIO(void): MeshIO()
{
}

MNIObjIO::MNIObjIO(const char *filename): MeshIO()
{
	read(filename);
}

MNIObjIO::~MNIObjIO(void)
{
}

void MNIObjIO::read(const char *filename)
{
	int nVertex = 0;
	int nFace = 0;
	int nNormal = 0;
	char buf[255];
	
	ifstream fin(filename);
	fin.getline(buf, sizeof(buf));
	if (buf[0] != 'P')
	{
		cout << "Unsupported input!\n";
		return;
	}

	sscanf(buf, "P %*f %*f %*f %*f %*f %d", &nVertex);
	nNormal = nVertex;
	for (int i = 0; i < nVertex * 2; )
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0) i++;
	}
	fin.getline(buf, sizeof(buf));
	fin.getline(buf, sizeof(buf));
	nFace = atoi(buf);

	if (nNormal != nVertex) nNormal = nVertex;
	else m_hasNormal = true;
	
	m_nVertex = nVertex;
	m_nFace = nFace;
	m_nNormal = nNormal;

	initArray();

	nVertex = 0;
	nFace = 0;
	nNormal = 0;

	fin.clear();
	fin.seekg(0, ios::beg);
	fin.getline(buf, sizeof(buf));
	for (; nVertex < m_nVertex; nVertex++)
	{
		fin.getline(buf, sizeof(buf));
		sscanf(buf, "%f %f %f", &m_vertex[nVertex * 3], &m_vertex[nVertex * 3 + 1], &m_vertex[nVertex * 3 + 2]);
	}
	fin.getline(buf, sizeof(buf));
	for (; nNormal < m_nNormal; nNormal++)
	{
		fin.getline(buf, sizeof(buf));
		sscanf(buf, "%f %f %f", &m_normal[nNormal * 3], &m_normal[nNormal * 3 + 1], &m_normal[nNormal * 3 + 2]);
	}
	fin.getline(buf, sizeof(buf));
	fin.getline(buf, sizeof(buf));
	fin.getline(buf, sizeof(buf));
	fin.getline(buf, sizeof(buf));
	fin.getline(buf, sizeof(buf));
	while (strlen(buf) > 0) fin.getline(buf, sizeof(buf));

	while (nFace < m_nFace * 6)
	{
		fin.getline(buf, sizeof(buf));
		char *tokens;
		char *ptr = strtok(buf, " ");
		while (ptr != NULL)
		{
			m_face[nFace++] = atoi(ptr);
			if (nFace % 3 == 0)
			{
				m_face[nFace] = m_face[nFace - 3]; nFace++;
				m_face[nFace] = m_face[nFace - 3]; nFace++;
				m_face[nFace] = m_face[nFace - 3]; nFace++;
			}
			ptr = strtok(NULL, " \r\n");
		}
	}
	
	fin.close();
}

