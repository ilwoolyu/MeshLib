/*************************************************
*	SurfaceUtil.h
*
*	Release: February 2015
*	Update: Jun 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>
#include "Mesh.h"

#ifndef SURFUTIL_HH_
#define SURFUTIL_HH_

class SurfaceUtil
{
public:
	// edge element
	struct edge
	{
		int vid1;
		int vid2;
		int fid1;
		int fid2;
		Vector tan;	// edge direction
		double len;	// edge length
		double T[3][3];	// tensor
		bool operator== (const edge& e) const
		{
			return vid1 == e.vid1 && vid2 == e.vid2;
		}
		bool operator< (const edge& e) const
		{
			if (vid1 != e.vid1)
				return vid1 < e.vid1;
			else
				return vid2 < e.vid2;
		}
	};

public:
	static void smoothing(Mesh *mesh, int iter);
	static void smoothing(Mesh *mesh, float sigma, int iter);
	static void smoothing(const Mesh *mesh, int iter, float *data);
	static void smoothing(const Mesh *mesh, int iter, double *data);
	static void curvature(const Mesh *mesh, double *cmin, double *cmax, double **umin = NULL, double **umax = NULL, int nIterSmoothingTensor = 0);
	static void curvature(const Mesh *mesh, float *cmin, float *cmax, float **umin = NULL, float **umax = NULL, int nIterSmoothingTensor = 0);
	static void sphere(Mesh *mesh, int type = 0, int max_iter = 10000);

private:
	static void smoothing(Mesh *mesh);
	static void smoothing(Mesh *mesh, float sigma);
	static void smoothing(const Mesh *mesh, float *data);
	static void smoothing(const Mesh *mesh, double *data);
	static void tensor(const Mesh *mesh, std::vector<double **> &Tv);
	static void smoothingTensor(const Mesh *mesh, const std::vector<double **> Tv, int iter);
	static void principalCurvature(const std::vector<double **> Tv, const int n, double *cmin, double *cmax, double **umin, double **umax);
};

#endif
