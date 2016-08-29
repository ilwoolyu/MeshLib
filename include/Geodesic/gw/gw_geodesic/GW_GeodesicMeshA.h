
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_GeodesicMesh.h
 *  \brief  Definition of class \c GW_GeodesicMesh
 *  \author Gabriel Peyré
 *  \date   4-9-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_GEODESICMESHA_H_
#define _GW_GEODESICMESHA_H_

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"
#include "../gw_core/GW_Face.h"
#include "../gw_core/GW_FaceIterator.h"
#include "../gw_core/GW_VertexIterator.h"
#include "GW_GeodesicVertex.h"
#include "GW_GeodesicFace.h"

namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_GeodesicMesh
 *  \brief  A mesh designed for computation of geodesic.
 *  \author Gabriel Peyré
 *  \date   4-9-2003
 *
 *  Overload class factory method to create \c GW_GeodesicVertex.
 */ 
/*------------------------------------------------------------------------------*/

class GW_GeodesicMeshA: public GW_Mesh
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    GW_GeodesicMeshA();
    virtual ~GW_GeodesicMeshA();
    //@}

	//-------------------------------------------------------------------------
	/** \name Class factory methods. */
	//-------------------------------------------------------------------------
	//@{
	virtual GW_Vertex& CreateNewVertex();
	virtual GW_Face& CreateNewFace();
	//@}

    //-------------------------------------------------------------------------
    /** \name Fast marching computations. */
    //-------------------------------------------------------------------------
	//@{
	void ResetGeodesicMesh();
	void ResetParametrizationData();
	void AddStartVertex( GW_GeodesicVertex& StartVert );
	void PerformFastMarching( GW_GeodesicVertex* pStartVertex=NULL );
	void SetUpFastMarching( GW_GeodesicVertex* pStartVertex=NULL );
	GW_Bool PerformFastMarchingOneStep();
	void PerformFastMarchingFlush();
	GW_Bool IsFastMarchingFinished();
    //@}

	void SetUseUnfolding( GW_Bool bUseUnfolding );
	GW_Bool GetUseUnfolding( );

    //-------------------------------------------------------------------------
    /** \name Callback management. */
    //-------------------------------------------------------------------------
    //@{
    typedef const double * (*T_V1CallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterV1CallbackFunction( T_V1CallbackFunction pFunc );
    typedef const double * (*T_V2CallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterV2CallbackFunction( T_V2CallbackFunction pFunc );
    typedef double (*T_Lam1CallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterLam1CallbackFunction( T_Lam1CallbackFunction pFunc );
    typedef double (*T_Lam2CallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterLam2CallbackFunction( T_Lam2CallbackFunction pFunc );
	
	typedef GW_Float (*T_WeightCallbackFunction)( GW_GeodesicVertex& Vert1, GW_Vector3D& Vert2 );
	void RegisterWeightCallbackFunction( T_WeightCallbackFunction pFunc );
	typedef GW_Bool (*T_FastMarchingCallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterForceStopCallbackFunction( T_FastMarchingCallbackFunction pFunc );
	typedef void (*T_NewDeadVertexCallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterNewDeadVertexCallbackFunction( T_NewDeadVertexCallbackFunction pFunc );
	typedef GW_Bool (*T_VertexInsersionCallbackFunction)( GW_GeodesicVertex& Vert, GW_Float rNewDist );
	void RegisterVertexInsersionCallbackFunction( T_VertexInsersionCallbackFunction pFunc );
	typedef GW_Float (*T_HeuristicToGoalCallbackFunction)( GW_GeodesicVertex& Vert );
	void RegisterHeuristicToGoalCallbackFunction( T_HeuristicToGoalCallbackFunction pFunc );
	//@}

	virtual GW_Vertex* GetRandomVertex( GW_Bool bForceFar = GW_True );

	static GW_Float BasicWeightCallback(GW_GeodesicVertex& Vert1, GW_Vector3D& Vert2);


protected:

	/** should be filled with the starting point of the marching before
	    calling PerformFastMarching */
	T_GeodesicVertexVector ActiveVertex_;

	/** a function that specify the metric tensor on the mesh */
	T_V1CallbackFunction V1Callback_;
	T_Lam1CallbackFunction Lam1Callback_;
	T_Lam2CallbackFunction Lam2Callback_;
	
	/** a function that specify the metric on the mesh */
	T_WeightCallbackFunction WeightCallback_;
	/** the callback function used to test if we should terminate the fast marching or not */
	T_FastMarchingCallbackFunction ForceStopCallback_;
	/** the callback function used to test when a new dead vertex is created */
	T_NewDeadVertexCallbackFunction NewDeadVertexCallback_;
	/** a function called to know if we should insert this vertex */
	T_VertexInsersionCallbackFunction VertexInsersionCallback_;
	/** a function called to give an heuristic for the remaining distance */
	T_HeuristicToGoalCallbackFunction HeuristicToGoalCallbackFunction_;
	
	/** just to controle interactive mode */
	GW_Bool bIsMarchingBegin_;
	GW_Bool bIsMarchingEnd_;


private:

	GW_Float ComputeVertexDistance( GW_GeodesicFace& CurrentFace, GW_GeodesicVertex& CurrentVertex, 
									GW_GeodesicVertex& Vert1, GW_GeodesicVertex& Vert2, GW_GeodesicVertex& CurrentFront );

	static GW_GeodesicVertex* UnfoldTriangle( GW_GeodesicFace& CurFace, GW_GeodesicVertex& v, GW_GeodesicVertex& v1, GW_GeodesicVertex& v2, GW_Float& dist, GW_Float& dot1, GW_Float& dot2);

	GW_Float ComputeUpdate_FastSweepingMethod( GW_GeodesicVertex& CurrentVertex, GW_GeodesicVertex& Vert1, GW_GeodesicVertex& Vert2, GW_Float d1, GW_Float d2, GW_Float dot, GW_Float F1, GW_Float F2 );
	GW_Float ComputeUpdate_SamplingMethod( GW_GeodesicVertex& CurrentVertex, GW_GeodesicVertex& Vert1, GW_GeodesicVertex& Vert2, GW_Float d1, GW_Float d2 );

	static GW_Float ComputeUpdate_SethianMethod( GW_Float d1, GW_Float d2, GW_Float a, GW_Float b, GW_Float dot, GW_Float F );
	static GW_Float ComputeUpdate_MatrixMethod( GW_Float d1, GW_Float d2, GW_Float a, GW_Float b, GW_Float dot, GW_Float F );

	/** Do we use unfolding to correct problem with non acute angles ? */
	static GW_Bool bUseUnfolding_;

};

/*------------------------------------------------------------------------------*/
/** \name a vector of GW_GeodesicMesh */
/*------------------------------------------------------------------------------*/
//@{
typedef std::vector<class GW_GeodesicMesh*> T_GeodesicMeshVector;
typedef T_GeodesicMeshVector::iterator IT_GeodesicMeshVector;
typedef T_GeodesicMeshVector::reverse_iterator RIT_GeodesicMeshVector;
typedef T_GeodesicMeshVector::const_iterator CIT_GeodesicMeshVector;
typedef T_GeodesicMeshVector::const_reverse_iterator CRIT_GeodesicMeshVector;
//@}


} // End namespace GW

#ifdef GW_USE_INLINE
    #include "GW_GeodesicMeshA.inl"
#endif


#endif // _GW_GEODESICMESH_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
