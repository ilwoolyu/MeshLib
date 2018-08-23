/*=================================================================
% perform_front_propagation_mesh - perform a Fast Marching front propagation on a 3D mesh.
%
%   [D,S,Q] = perform_front_propagation_mesh(vertex, faces, W,start_points,end_points, nb_iter_max,H,L, values, dmax);
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 2 x num_start_points matrix where k is the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

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

#include <math.h>
#include <algorithm>
#include <map>
#include <list>
#include <string>
#include "GeodesicA.h"
#include "SurfaceUtil.h"
#include "Util/Geom.h"

GeodesicA::GeodesicA(const Mesh *mesh)
{
	nverts = -1;
	nfaces = -1;
	
	V1 = NULL;	// principal direction 1 of the metric tensor
	V1_ = NULL;	// workspace for principal direction 1 of the metric tensor
	Lam1 = NULL;	// eigenvalue 1 of the metric tensor
	Lam2 = NULL;	// eigenvalue 2 of the metric tensor

	D = NULL;	// distance
	S = NULL;	// state
	Q = NULL;	// nearest neighbor
	A = NULL;	// cumulative area
	T = NULL;	// tensor
	pT = NULL;	// work space for tensor
	exc	= NULL;
	varea = NULL;
	
	if (mesh != NULL) setupMesh(mesh);
	
	// for callback functions
	setInstance(this);
}

GeodesicA::~GeodesicA(void)
{
	if (V1 != NULL) delete [] V1;
	if (V1_ != NULL) delete [] V1_;
	if (Lam1 != NULL) delete [] Lam1;
	if (Lam2 != NULL) delete [] Lam2;
	if (varea != NULL) delete [] varea;
	
	if (D != NULL) delete [] D;
	// second output : state
	if (S != NULL) delete [] S;
	// second output : segmentation
	if (Q != NULL) delete [] Q;
	if (A != NULL) delete [] A;
}

GeodesicA *GeodesicA::instance = 0;
void GeodesicA::setInstance(GeodesicA *_instance)
{
	instance = _instance;
}

const double * GeodesicA::V1Callback(GW_GeodesicVertex& Vert1)
{
	GW_U32 i = Vert1.GetID();
	return instance->V1[i];
}

double GeodesicA::Lam1Callback(GW_GeodesicVertex& Vert1)
{
	GW_U32 i = Vert1.GetID();
	return instance->Lam1[i];
}

double GeodesicA::Lam2Callback(GW_GeodesicVertex& Vert1)
{
	GW_U32 i = Vert1.GetID();
	return instance->Lam2[i];
}

GW_Float GeodesicA::WeightCallback(GW_GeodesicVertex& Vert1, GW_Vector3D& Vert2)
{
	GW_U32 i = Vert1.GetID();
	
	GW_Vector3D E = Vert2 - Vert1.GetPosition();
	GW_Vector3D N = Vert1.GetNormal();
	GW_Float p = N * E;
	E -= N * p;
	E /= E.Norm();
	
	const double **T = (const double **) instance->T;
	GW_Float speed = E[0] * (T[i][0] * E[0] + T[i][1] * E[1] + T[i][2] * E[2]) +
					E[1] * (T[i][3] * E[0] + T[i][4] * E[1] + T[i][5] * E[2]) +
					E[2] * (T[i][6] * E[0] + T[i][7] * E[1] + T[i][8] * E[2]);

	/*if (false)
	{
		cout << i << endl;
		cout << Vert1.GetPosition()[0] << " " << Vert1.GetPosition()[1] << " " << Vert1.GetPosition()[2] << ": ";
		cout << Vert2[0] << " " << Vert2[1] << " " << Vert2[2] << ": ";
		cout << speed << endl;
		cout << T[i][0] << " " << T[i][1] << " " << T[i][2] << endl;
		cout << T[i][3] << " " << T[i][4] << " " << T[i][5] << endl;
		cout << T[i][6] << " " << T[i][7] << " " << T[i][8] << endl;
		cout << endl;
	}*/
//	return instance->Ww[i] / speed;
	return 1 / speed;
}

GW_Bool GeodesicA::StopMarchingCallback( GW_GeodesicVertex& Vert )
{
	// check if the end point has been reached
	GW_U32 i = Vert.GetID();
	//printf("ind %d\n",i );
	//printf("dist %f\n",Vert.GetDistance() );
	if( Vert.GetDistance() > instance->dmax )
		return true;
	if (instance->nbr_iter >= instance->niter_max)
		return true;
	if (instance->area_max > 0)
	{
		instance->area_prop += instance->varea[i];
		instance->A[i] = instance->area_prop;
		if (instance->area_prop > instance->area_max)
			return true;
	}
	
	for( int k=0; k<instance->nend; ++k )
		if( instance->end_points[k]==i )
			return true;
	return false;
}

GW_Bool GeodesicA::InsersionCallback( GW_GeodesicVertex& Vert, GW_Float rNewDist )
{
	// check if the distance of the new point is less than the given distance
	GW_U32 i = Vert.GetID();
	bool doinsersion = instance->nbr_iter <= instance->niter_max;
	if( instance->L!=NULL )
		doinsersion = doinsersion && (rNewDist<instance->L[i]);
	if ( instance->exc != NULL)
		doinsersion = doinsersion && !instance->exc[i];
	instance->nbr_iter++;
	return doinsersion;
}

GW_Float GeodesicA::HeuristicCallback( GW_GeodesicVertex& Vert )
{
	// return the heuristic distance
	GW_U32 i = Vert.GetID();
	return instance->H[i];
}

void GeodesicA::setupCurvature(const Mesh *mesh)
{
	int nv = mesh->nVertex();
	m_cmin = new double[nv];
	m_cmax = new double[nv];
	m_umin = new double*[nv]; for (int i = 0; i < nv; i++) m_umin[i] = new double[3];
	m_umax = new double*[nv]; for (int i = 0; i < nv; i++) m_umax[i] = new double[3];
	SurfaceUtil::curvature(mesh, m_cmin, m_cmax, m_umin, m_umax);
}

void GeodesicA::setupMesh(const Mesh *mesh)
{
	setupCurvature(mesh);
	
	varea = new double[mesh->nVertex()];
	for (int i = 0; i < mesh->nVertex(); i++)
		varea[i] = vertexArea(mesh, i);
	// arg1 : vertex
	nverts = mesh->nVertex();
	// arg2 : faces
	nfaces = mesh->nFace();

	// create the mesh
	GWMesh.SetNbrVertex(nverts);
	for(int i = 0; i < nverts; ++i)
	{
		GW_GeodesicVertex& vert = (GW_GeodesicVertex&) GWMesh.CreateNewVertex();
		vert.SetPosition(GW_Vector3D((*mesh->vertex(i))[0],(*mesh->vertex(i))[1],(*mesh->vertex(i))[2]));
		GWMesh.SetVertex(i, &vert);
	}
	GWMesh.SetNbrFace(nfaces);
	for( int i = 0; i < nfaces; ++i)
	{
		GW_GeodesicFace& face = (GW_GeodesicFace&) GWMesh.CreateNewFace();
		GW_Vertex* v1 = GWMesh.GetVertex((*mesh->face(i))[0]); GW_ASSERT(v1 != NULL);
		GW_Vertex* v2 = GWMesh.GetVertex((*mesh->face(i))[1]); GW_ASSERT(v2 != NULL);
		GW_Vertex* v3 = GWMesh.GetVertex((*mesh->face(i))[2]); GW_ASSERT(v3 != NULL);
		face.SetVertex( *v1,*v2,*v3 );
		GWMesh.SetFace(i, &face);
	}
	GWMesh.BuildConnectivity();
	
	for(int i = 0; i < nverts; ++i)
	{
		GW_Vertex* v = GWMesh.GetVertex(i);
		v->BuildRawNormal();
	}
	
	T = new double*[nverts];
	pT = new double[nverts * 9];
	for (int i = 0; i < nverts; i++)
	{
		T[i] = &pT[i * 9];
		T[i] = &pT[i * 9];
	}

	if (V1 != NULL) delete [] V1;
	V1 = new double*[nverts];
	if (V1_ != NULL) delete [] V1_;
	V1_ = new double[nverts * 3];
	if (Lam1 != NULL) delete [] Lam1;
	Lam1 = new double[nverts];
	if (Lam2 != NULL) delete [] Lam2;
	Lam2 = new double[nverts];
	for (int i = 0; i < nverts; i++) V1[i] = &V1_[i * 3];

	GWMesh.RegisterV1CallbackFunction(V1Callback);
	GWMesh.RegisterLam1CallbackFunction(Lam1Callback);
	GWMesh.RegisterLam2CallbackFunction(Lam2Callback);
	
	GWMesh.RegisterWeightCallbackFunction(WeightCallback);
	GWMesh.RegisterForceStopCallbackFunction(StopMarchingCallback);
	GWMesh.RegisterVertexInsersionCallbackFunction(InsersionCallback);
	
	setupCurvatureTensor();

	// first ouput : distance
	if (D != NULL) delete [] D;
	D = new double[nverts];
	// second output : state
	if (S != NULL) delete [] S;
	S = new int[nverts];
	if (Q != NULL) delete [] Q;
	// second output : segmentation, nearest neighbor
	Q = new int[nverts];
	if (A != NULL) delete [] A;
	// second output : cumulative area
	A = new double[nverts];

	setupOptions();
}

void GeodesicA::setArea(const double *area)
{
	for (int i = 0; i < nverts; i++)
		varea[i] = area[i];
}

void GeodesicA::setExclusion(const bool *list)
{
	if (exc == NULL) exc = new bool[nverts];
	memcpy(exc, list, sizeof(bool) * nverts);
}

void GeodesicA::setTensor(const double **_T)
{
	for (int i = 0; i < nverts; i++)
	{
		memcpy(T[i], _T[i], sizeof(double) * 9);
	}

	// eigendecomposition of tensor metrics
	// warning: make sure your tensor should be rank 2 matrix
	for (int i = 0; i < nverts; i++)
	{
		double T1[3][3] = {{T[i][0], T[i][1], T[i][2]}, {T[i][3], T[i][4], T[i][5]}, {T[i][6], T[i][7], T[i][8]}};
	
		double lambda[3];
		double eigv[3][3];
		
		LinearAlgebra::eig3symmetric(T1, lambda, eigv);
		//GW_ASSERT(fabs(lambda[0]) < GW_EPSILON); 	// lambda0 should be always zero
		//cout << "lambda[0]: " << lambda[0] << endl;
		V1[i][0] = eigv[2][0]; V1[i][1] = eigv[2][1]; V1[i][2] = eigv[2][2]; // the first eigenvector
		Lam1[i] = lambda[2];
		Lam2[i] = lambda[1];
	}
}

void GeodesicA::setTensor(const double **_u, const double *_lam1, const double *_lam2)
{
	GW_Vector3D *V2 = new GW_Vector3D[nverts];
	for (int i = 0; i < nverts; i++)
	{
		GW_Vector3D V1_(_u[i][0], _u[i][1], _u[i][2]);
		GW_Vector3D N = GWMesh.GetVertex(i)->GetNormal();
		
		// force orthogonal to normal
		double p = N * V1_;
		V1_ -= N * p; V1_ /= V1_.Norm();
		
		V2[i] = N ^ V1_; V2[i] /= V2[i].Norm();
		V1[i][0] = V1_[0]; V1[i][1] = V1_[1]; V1[i][2] = V1_[2];
		
		Lam1[i] = _lam1[i];
		Lam2[i] = _lam2[i];
	}

	for (int i = 0; i < nverts; i++)
	{
		// covariance
		for (int j = 0; j < 3; j++)
			for (int k = j; k < 3; k++)
				T[i][j * 3 + k] = T[i][k * 3 + j] = Lam1[i] * V1[i][j] * V1[i][k] + Lam2[i] * V2[i][j] * V2[i][k];
		/*if (i == 194)
		{
			cout << Lam1[i] << " " << Lam2[i] << endl;
			cout << "u: " << _u[i][0] << " " << _u[i][1] << " " << _u[i][2] << endl;
			cout << "V1: " << V1[i][0] << " " << V1[i][1] << " " << V1[i][2] << endl;
			cout << "V2: " << V2[i] << endl << endl;
			cout << T[i][0] << " " << T[i][1] << " " << T[i][2] << endl;
			cout << T[i][3] << " " << T[i][4] << " " << T[i][5] << endl;
			cout << T[i][6] << " " << T[i][7] << " " << T[i][8] << endl<<endl;;
		}*/
	}
	
	delete [] V2;
}

void GeodesicA::setupCurvatureTensor(void)
{
	double **_u = new double*[nverts];
	double *_pu = new double[nverts * 3];
	for (int i = 0; i < nverts; i++) _u[i] = &_pu[i * 3];
	double *_lam1 = new double[nverts];
	double *_lam2 = new double[nverts];
	
	double eps = 1e-2;
	double m = FLT_MAX, M = 0;
	for (int i = 0; i < nverts; i++)
	{
		GW_Vertex *v = GWMesh.GetVertex(i);
		//v->BuildCurvatureData();
		//GW_Float cmin = v->GetMinCurv(), cmax = v->GetMaxCurv();
		GW_Float cmin = m_cmin[i], cmax = m_cmax[i];
		
		GW_Float invcmin = 1 / (fabs(cmin) + eps);
		GW_Float invcmax = 1 / (fabs(cmax) + eps);
		if (m > invcmin) m = invcmin;
		if (M < invcmin) M = invcmin;
		if (m > invcmax) m = invcmax;
		if (M < invcmax) M = invcmax;
	}
	
	cout << "Min Curv: " << m << " ";
	cout << "Max Curv: " << M << endl;
	GW_Float a = 1, b = 100;
	for (int i = 0; i < nverts; i++)
	{
		GW_Vertex *v = GWMesh.GetVertex(i);
		//GW_Vector3D dmin = v->GetMinCurvDirection(), dmax = v->GetMaxCurvDirection();
		//GW_Float cmin = v->GetMinCurv(), cmax = v->GetMaxCurv();
		const double *dmin = m_umin[i], *dmax = m_umax[i];
		GW_Float cmin = m_cmin[i], cmax = m_cmax[i];
		
		GW_Float invcmin = 1 / (fabs(cmin) + eps);
		GW_Float invcmax = 1 / (fabs(cmax) + eps);
		GW_Float recmin = (b - a) * (invcmin - m) / (M - m) + a;
		GW_Float recmax = (b - a) * (invcmax - m) / (M - m) + a;
		_lam1[i] = recmin;
		_lam2[i] = recmax;
		_u[i][0] = dmin[0]; _u[i][1] = dmin[1]; _u[i][2] = dmin[2];
	}
	setTensor((const double **)_u, (const double *)_lam1, (const double *)_lam2);
	
	/*FILE *fp = fopen("/NIRAL/work/ilwoolyu/Surfaces/Gyrification/test/eig.txt", "w");
	for (int i = 0; i < nverts; i++)
	{
		double lambda[3];
		double eigv[3][3];
		double T1[3][3] = {{T[i][0], T[i][1], T[i][2]}, {T[i][3], T[i][4], T[i][5]}, {T[i][6], T[i][7], T[i][8]}};
		LinearAlgebra::eig3symmetric(T1, lambda, eigv);
		//fprintf(fp, "%f %f %f\n", eigv[2][0], eigv[2][1], eigv[2][2]);
		fprintf(fp, "%f %f %f\n", _u[i][0], _u[i][1], _u[i][2]);
	}
	fclose(fp);*/
	
	delete [] _u;
	delete [] _pu;
	delete [] _lam1;
	delete [] _lam2;
}

double GeodesicA::vertexArea(const Mesh *mesh, int id)
{
	const int *list = mesh->vertex(id)->list();
	int nn = mesh->vertex(id)->nNeighbor();

	double sum = 0;
	for (int j = 0; j < nn; j++)
	{
		int id1 = mesh->vertex(id)->list(j);
		int id2 = mesh->vertex(id)->list((j + 1) % nn);
		Vector E1 = Vector(mesh->vertex(id)->fv(), mesh->vertex(id1)->fv());
		Vector E2 = Vector(mesh->vertex(id)->fv(), mesh->vertex(id2)->fv());
		double a = E1.norm(), b = E2.norm();
		E1.unit(); E2.unit();
		double inner = E1 * E2;
		if (inner > 1) inner = 1;
		if (inner < -1) inner = -1;
		double theta = acos(inner);
		double area = a * b * sin(theta) / 2;
		sum += area / 3;
	}
	
	return sum;
}

void GeodesicA::setupOptions(const double *_Ww, const double *_H, const double *_L, const double *_values)
{
	// arg3 : Ww - weight
	Ww = _Ww;
	if (Ww == NULL)
	{
		double *W = new double[nverts];
		for (int i = 0; i < nverts; i++)
			W[i] = 1;
		Ww = (const double *)W;
	}

	// arg7 : H
	H = _H;
	// arg8 : L
	L = _L;
	// argument 9: value list
	values = _values;

	if (H != NULL)
		GWMesh.RegisterHeuristicToGoalCallbackFunction(HeuristicCallback);
	else
		GWMesh.RegisterHeuristicToGoalCallbackFunction(NULL);
	// initialize the distance of the starting points
	if (values != NULL)
	{
		for (int i = 0; i < nstart; ++i)
		{
			GW_GeodesicVertex* v = (GW_GeodesicVertex*) GWMesh.GetVertex((GW_U32) start_points[i]);
			GW_ASSERT(v != NULL);
			v->SetDistance(values[i]);
		}
	}
}

void GeodesicA::perform_front_propagation(int start_point, int end_point)
{
	int start_points[1] = {start_point};
	dmax = 1e9;
	
	if (end_point == -1)
	{
		perform_front_propagation(start_points, 1, NULL, 0);
	}
	else
	{
		int end_points[1] = {end_point};
		perform_front_propagation(start_points, 1, end_points, 1);
	}
}

void GeodesicA::perform_front_propagation(int start_point, double _dmax)
{
	int start_points[1] = {start_point};
	
	perform_front_propagation(start_points, 1, NULL, 0, _dmax);
}

void GeodesicA::perform_front_propagation(const int *_start_points, int _nstart, const int *_end_points, int _nend, double _dmax, int _niter_max, double _area_max)
{
	nbr_iter = 0;
	area_prop = 0;

	area_max = _area_max;
	nstart = _nstart;
	nend = _nend;
	
	// arg4 : start_points
	start_points = _start_points;
	// arg5 : end_points
	end_points = _end_points;
	// arg6 : niter_max
	niter_max = (_niter_max == 0 || nverts * 1.2 < _niter_max)? nverts * 1.2: _niter_max;

	// argument 10: dmax
	dmax = _dmax;

	// set up fast marching
	GWMesh.ResetGeodesicMesh();
	for(int i = 0; i < nstart; ++i)
	{
		GW_GeodesicVertex* v = (GW_GeodesicVertex*)GWMesh.GetVertex((GW_U32)start_points[i]);
		GW_ASSERT(v != NULL);
		GWMesh.AddStartVertex(*v);
	}
	
	GWMesh.SetUpFastMarching();
	
	// perform fast marching
	GWMesh.PerformFastMarching();
	
	// output result
	for(int i = 0; i < nverts; ++i)
	{
		GW_GeodesicVertex* v = (GW_GeodesicVertex*)GWMesh.GetVertex((GW_U32)i);
		GW_ASSERT(v != NULL);
		D[i] = v->GetDistance();
		S[i] = v->GetState();
		GW_GeodesicVertex* v1 = v->GetFront();
		if(v1 == NULL)
			Q[i] = -1;
		else
			Q[i] = v1->GetID();
	}
}

const double * GeodesicA::dist(void)
{
	return D;
}

const int * GeodesicA::state(void)
{
	return S;
}

const int * GeodesicA::source(void)
{
	return Q;
}

const double * GeodesicA::area(void)
{
	return A;
}

int GeodesicA::niter(void)
{
	return nbr_iter;
}
