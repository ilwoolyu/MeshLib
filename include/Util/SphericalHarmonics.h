/*************************************************
*	SphericalHarmonics.h
*
*	Release: August 2012
*	Update: June 2021
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*	
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
*************************************************/

#pragma once
#include "Geom.h"

class SphericalHarmonics
{
public:
	static void basis(int degree, float *p, float *Y, int from = 0);
	static void basis(int degree, double *p, double *Y, int from = 0);
};
