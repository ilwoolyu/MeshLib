/*************************************************
*	Geodesic distance wrapper for MeshLib C++
*
*	Release: April 2015
*	Update: July 2016
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>
#include "Mesh.h"

#ifndef GEODESICA_HH_
#define GEODESICA_HH_

using std::string;

#include "gw/config.h"
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "gw/gw_geodesic/GW_GeodesicMeshA.h"

using namespace GW;

class GeodesicA
{
public:
	GeodesicA(const Mesh *mesh = NULL);
	~GeodesicA(void);
	void setupMesh(const Mesh *mesh);
	void setupOptions(const double *Ww = NULL, const double *H = NULL, const double *L = NULL, const double *values = NULL);
	void perform_front_propagation(int start_point, int end_point = -1);
	void perform_front_propagation(int start_point, double dmax);
	void perform_front_propagation(const int *start_points, int nstart, const int *end_points, int nend, double dmax = 1e9, int niter_max = 0, double _area_max = 0);
	void setTensor(const double **_T);
	void setTensor(const double **_u, const double *_lam1, const double *_lam2);
	void setExclusion(const bool *list);
	void setArea(const double *area);
	const double *dist(void);
	const int *state(void);
	const int *source(void);
	const double *area(void);
	int niter(void);
	
private:
	static void setInstance(GeodesicA *_instance);
	
	static const double * V1Callback(GW_GeodesicVertex& Vert1);
	static double Lam1Callback(GW_GeodesicVertex& Vert1);
	static double Lam2Callback(GW_GeodesicVertex& Vert1);
	
	static GW_Float WeightCallback(GW_GeodesicVertex& Vert1, GW_Vector3D& Vert2);
	static GW_Bool StopMarchingCallback(GW_GeodesicVertex& Vert);
	static GW_Bool InsersionCallback(GW_GeodesicVertex& Vert, GW_Float rNewDist);
	static GW_Float HeuristicCallback(GW_GeodesicVertex& Vert);
	void setupCurvature(const Mesh *mesh);
	void setupCurvatureTensor(void);
	double vertexArea(const Mesh *mesh, int id);

private:
	int nverts; 
	int nfaces;
	int nstart, nend;
	int niter_max;
	double dmax;	
	const int *start_points, *end_points;
	const double *Ww, *H, *L, *values;
	double **T, *pT;
	// curvature
	double *m_cmin;
	double *m_cmax;
	double **m_umin;
	double **m_umax;
	
	// metric tensor
	double **V1;
	double *V1_;
	double *Lam1, *Lam2;
	
	bool *exc;	// excluded points
	double *varea;	// vertex area

	// outputs 
	double *D;	// distance
	int *S;	// state
	int *Q;	// nearest neighbor
	double *A;	// cumulative area
	int nbr_iter;
	double area_max;	// max area
	double area_prop;	// propagation area
	static GeodesicA *instance;
	GW_GeodesicMeshA GWMesh;
};

#endif
