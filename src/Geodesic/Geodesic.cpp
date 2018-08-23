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
*	Update: February 2016
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

#include "Geodesic.h"

Geodesic::Geodesic(const Mesh *mesh)
{
	nverts = -1;
	nfaces = -1;

	D = NULL;	// distance
	S = NULL;	// state
	Q = NULL;	// nearest neighbor
	
	if (mesh != NULL) setupMesh(mesh);
	
	// for callback functions
	setInstance(this);
}

Geodesic::~Geodesic(void)
{
	if (D != NULL) delete [] D;
	// second output : state
	if (S != NULL) delete [] S;
	// second output : segmentation
	if (Q != NULL) delete [] Q;
}

Geodesic *Geodesic::instance = 0;
void Geodesic::setInstance(Geodesic *_instance)
{
	instance = _instance;
}

GW_Float Geodesic::WeightCallback(GW_GeodesicVertex& Vert)
{
	GW_U32 i = Vert.GetID();
	return instance->Ww[i];
}

GW_Bool Geodesic::StopMarchingCallback( GW_GeodesicVertex& Vert )
{
	// check if the end point has been reached
	GW_U32 i = Vert.GetID();
//	display_message("ind %d",i );
//	display_message("dist %f",Vert.GetDistance() );
	if( Vert.GetDistance() > instance->dmax )
		return true;
	for( int k=0; k<instance->nend; ++k )
		if( instance->end_points[k]==i )
			return true;
	return false;
}

GW_Bool Geodesic::InsersionCallback( GW_GeodesicVertex& Vert, GW_Float rNewDist )
{
	// check if the distance of the new point is less than the given distance
	GW_U32 i = Vert.GetID();
	bool doinsersion = instance->nbr_iter <= instance->niter_max;
	if( instance->L!=NULL )
		doinsersion = doinsersion && (rNewDist<instance->L[i]);
	instance->nbr_iter++;
	return doinsersion;
}
GW_Float Geodesic::HeuristicCallback( GW_GeodesicVertex& Vert )
{
	// return the heuristic distance
	GW_U32 i = Vert.GetID();
	return instance->H[i];
}

void Geodesic::setupMesh(const Mesh *mesh)
{
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
	
	GWMesh.RegisterWeightCallbackFunction(WeightCallback);
	GWMesh.RegisterForceStopCallbackFunction(StopMarchingCallback);
	GWMesh.RegisterVertexInsersionCallbackFunction(InsersionCallback);
	
	// first ouput : distance
	if (D != NULL) delete [] D;
	D = new double[nverts];
	// second output : state
	if (S != NULL) delete [] S;
	S = new int[nverts];
	if (Q != NULL) delete [] Q;
	// second output : segmentation, nearest neighbor
	Q = new int[nverts];
	
	setupOptions();
}

void Geodesic::setupOptions(const double *_Ww, const double *_H, const double *_L, const double *_values)
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

void Geodesic::perform_front_propagation(int start_point, int end_point)
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

void Geodesic::perform_front_propagation(int start_point, double _dmax)
{
	int start_points[1] = {start_point};
	
	perform_front_propagation(start_points, 1, NULL, 0, _dmax);
}

void Geodesic::perform_front_propagation(const int *_start_points, int _nstart, const int *_end_points, int _nend, double _dmax, int _niter_max)
{
	nbr_iter = 0;

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

const double * Geodesic::dist(void)
{
	return D;
}

const int * Geodesic::state(void)
{
	return S;
}

const int * Geodesic::source(void)
{
	return Q;
}

