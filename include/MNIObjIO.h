/*************************************************
*	MNIObjIO.h
*
*	Release: March 2015
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#pragma once
#include "MeshIO.h"

using namespace std;

class MNIObjIO: public MeshIO
{
public:
	MNIObjIO(void);
	MNIObjIO(const char *filename);
	~MNIObjIO(void);
	void read(const char *filename);
};
