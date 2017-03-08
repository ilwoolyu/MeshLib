/*************************************************
*	Geodesic distance wrapper for MeshLib C++
*
*	Release: April 2015
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>

#ifndef GEODESIC_HH_
#define GEODESIC_HH_

using std::string;

#include "Mesh.h"
#include "gw/config.h"
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "gw/gw_geodesic/GW_GeodesicMesh.h"

using namespace GW;

class Geodesic
{
public:
	Geodesic(const Mesh *mesh = NULL);
	~Geodesic(void);
	void setupMesh(const Mesh *mesh);
	void setupOptions(const double *Ww = NULL, const double *H = NULL, const double *L = NULL, const double *values = NULL);
	void perform_front_propagation(int start_point, int end_point = -1);
	void perform_front_propagation(int start_point, double dmax);
	void perform_front_propagation(const int *start_points, int nstart, const int *end_points, int nend, double dmax = 1e9, int niter_max = 0);
	const double *dist(void);
	const int *state(void);
	const int *source(void);
	
private:
	static void setInstance(Geodesic *_instance);
	static GW_Float WeightCallback(GW_GeodesicVertex& Vert);
	static GW_Bool StopMarchingCallback(GW_GeodesicVertex& Vert);
	static GW_Bool InsersionCallback(GW_GeodesicVertex& Vert, GW_Float rNewDist);
	static GW_Float HeuristicCallback(GW_GeodesicVertex& Vert);

private:
	int nverts; 
	int nfaces;
	int nstart, nend;
	int niter_max;
	double dmax;	
	const int *start_points, *end_points;
	const double *Ww, *H, *L, *values;
	// outputs 
	double *D;	// distance
	int *S;	// state
	int *Q;	// nearest neighbor
	int nbr_iter;
	static Geodesic *instance;
	GW_GeodesicMesh GWMesh;
};

#endif
