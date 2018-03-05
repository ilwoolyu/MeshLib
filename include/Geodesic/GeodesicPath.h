/*=================================================================
%   path is a 3D curve that is the shortest path starting at x.
%   You can force to use a fully discrete descent using
%   options.method = 'discrete'.
%
%   Copyright(c) 2007 Gabriel Peyre
*=================================================================*/

/*************************************************
*	Geodesic path for MeshLib C++
*
*	Release: March 2018
*	Update: March 2018
*
*	Vanderbilt University
*	Electrical Engineering and Computer Science
*
*	Ilwoo Lyu, ilwoo.lyu@vanderbilt.edu
*************************************************/

#ifndef GEODESIC_PATH_HH_
#define GEODESIC_PATH_HH_

#include <stdio.h>
#include <vector>
#include "Mesh.h"

using namespace std;

class GeodesicPath
{
protected:
	class Edge
	{
	public:
		int vid1;
		int vid2;
		int face;

		Edge(void);
		Edge(int vid1, int vid2);
		Edge(int vid1, int vid2, int face);
		bool operator ==(const Edge &p) const;
		bool operator <(const Edge &p) const;
		bool operator >(const Edge &p) const;
	};
	vector<Edge> m_table;
	vector<Edge> m_vlist;
	vector<double> m_plist;
	vector<float *> m_path;
	const Mesh *m_mesh;
	const double *m_dist;

public:
	GeodesicPath(void);
	GeodesicPath(const double *dist, const Mesh *mesh);
	~GeodesicPath(void);
	const float *getPoint(int id);
	int size(void);
	void computeGeodesicPath(int x);
private:
	void init(void);
	void clear(void);
	int get_face_face(int f, int v, int w);
	int get_edge(int v, int w);
	int get_vertex_face(int f, int v, int w);
	void vertex_stepping(int v, int &w, int &f);
};
#endif
