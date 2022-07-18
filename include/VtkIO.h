/*************************************************
*	VtkIO.h
*
*	Release: November 2013
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#pragma once
#include <algorithm>
#include "MeshIO.h"

#if defined(_WIN32) || defined(_WIN64)
#define strcasecmp _stricmp
#endif

using namespace std;

class VtkIO: public MeshIO
{
public:
	VtkIO(void);
	VtkIO(const char *filename);
	~VtkIO(void);
	void read(const char *filename);
	static void save(const char *filename, const Mesh *mesh, bool normal = false, bool binary = false);
};
