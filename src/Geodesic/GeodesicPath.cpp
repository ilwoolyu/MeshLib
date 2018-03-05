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

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "Mesh.h"
#include "GeodesicPath.h"

using namespace std;

GeodesicPath::Edge::Edge(void)
{
	this->vid1 = 0;
	this->vid2 = 0;
	this->face = 0;
}
GeodesicPath::Edge::Edge(int vid1, int vid2)
{
	this->vid1 = vid1;
	this->vid2 = vid2;
	this->face = 0;
}
GeodesicPath::Edge::Edge(int vid1, int vid2, int face)
{
	this->vid1 = vid1;
	this->vid2 = vid2;
	this->face = face;
}
bool GeodesicPath::Edge::operator ==(const Edge &p) const
{
	return (this->vid1 == p.vid1) && (this->vid2 == p.vid2);
}
bool GeodesicPath::Edge::operator <(const Edge &p) const
{
	return (this->vid1 < p.vid1) || ((this->vid1 == p.vid1) && (this->vid2 < p.vid2));
}
bool GeodesicPath::Edge::operator >(const Edge &p) const
{
	return (this->vid1 > p.vid1) || ((this->vid1 == p.vid1) && (this->vid2 > p.vid2));
}
GeodesicPath::GeodesicPath(void)
{
	m_mesh = NULL;
}
GeodesicPath::GeodesicPath(const double *dist, const Mesh *mesh)
{
	m_mesh = mesh;
	m_dist = dist;
	init();
}
GeodesicPath::~GeodesicPath(void)
{
	m_mesh = NULL;
	clear();
}
const float * GeodesicPath::getPoint(int id)
{
	return m_path[id];
}
int GeodesicPath::size(void)
{
	return m_path.size();
}
void GeodesicPath::computeGeodesicPath(int x)
{
	clear();
	int w, f;
	vertex_stepping(x, w, f);
	m_vlist.push_back(m_table[get_edge(x, w)]);
	m_plist.push_back(1.0);
	float *p = new float[3];
	p[0] = m_mesh->vertex(x)->fv()[0];
	p[1] = m_mesh->vertex(x)->fv()[1];
	p[2] = m_mesh->vertex(x)->fv()[2];
	m_path.push_back(p);
	double Dprev = m_dist[x];

	while (true)
	{
		int n = m_vlist.size() - 1;
		int i = m_vlist[n].vid1;
		int j = m_vlist[n].vid2;
		int k = get_vertex_face(f, i, j);
		double a = m_dist[i], b = m_dist[j], c = m_dist[k];

		int f1 = get_face_face(f, i, k);
		int f2 = get_face_face(f, j, k);

		double x = m_plist[n];
		double y = 1 - x;

		Vector u = Vector(m_mesh->vertex(k)->fv(), m_mesh->vertex(i)->fv());
		Vector v = Vector(m_mesh->vertex(k)->fv(), m_mesh->vertex(j)->fv());
		
		double det = (u * u) * (v * v) - (u * v) * (u * v);
		double Ainv[2][2] = { {(v * v) / det, -(u * v) / det },{ -(u * v) / det, (u * u) / det } };
		double nx = Ainv[0][0] * (a - c) + Ainv[0][1] * (b - c);
		double ny = Ainv[1][0] * (a - c) + Ainv[1][1] * (b - c);

		double etas = -y / ny;
		double etat = -x / nx;
		double s = x + etas * nx;
		double t = y + etat * ny;

		if (etas < 0 && s >= 0 && s <= 1 && f1 >= 0)
		{
			m_plist.push_back(s);
			m_vlist.push_back(m_table[get_edge(i, k)]);
			f = f1;
		}
		else if (etat < 0 && t >= 0 && t <= 1 && f2 >= 0)
		{
			m_plist.push_back(t);
			m_vlist.push_back(m_table[get_edge(j, k)]);
			f = f2;
		}
		else
		{
			int z = (a <= b) ? i : j;
			vertex_stepping(z, w, f);
			m_vlist.push_back(m_table[get_edge(z, w)]);
			m_plist.push_back(1.0);
		}
		n++;
		double Dnew = m_dist[m_vlist[n].vid1] * m_plist[n] + m_dist[m_vlist[n].vid2] * (1 - m_plist[n]);
		Vector p = Vector(m_mesh->vertex(m_vlist[n].vid1)->fv()) * m_plist[n] + Vector(m_mesh->vertex(m_vlist[n].vid2)->fv()) * (1 - m_plist[n]);
		float *q = new float[3];
		q[0] = p.fv()[0];
		q[1] = p.fv()[1];
		q[2] = p.fv()[2];
		m_path.push_back(q);
		if (Dnew == 0 || (Dprev == Dnew && n > 0)) break;
		
		Dprev = Dnew;
	}
}
void GeodesicPath::init(void)
{
	for (int i = 0; i < m_mesh->nFace(); i++)
	{
		const int *id = m_mesh->face(i)->list();
		Edge e[3];
		e[0].vid1 = id[0]; e[0].vid2 = id[1]; e[0].face = i; m_table.push_back(e[0]);
		e[1].vid1 = id[1]; e[1].vid2 = id[2]; e[1].face = i; m_table.push_back(e[1]);
		e[2].vid1 = id[2]; e[2].vid2 = id[0]; e[2].face = i; m_table.push_back(e[2]);
	}
	sort(m_table.begin(), m_table.end());
}
void GeodesicPath::clear(void)
{
	m_vlist.clear();
	m_plist.clear();
	for (int i = 0; i < m_path.size(); i++)
		delete [] m_path[i];
	m_path.clear();
}
int GeodesicPath::get_face_face(int f, int v, int w)
{
	int f1 = m_table[get_edge(v, w)].face;
	int f2 = m_table[get_edge(w, v)].face;

	if (f == f1)
		return f2;
	else
		return f1;
}
int GeodesicPath::get_edge(int v, int w)
{
	int key = lower_bound(m_table.begin(), m_table.end(), Edge(v, w)) - m_table.begin();
	return key;
}
int GeodesicPath::get_vertex_face(int f, int v, int w)
{
	int id;
	const int *vlist = m_mesh->face(f)->list();
	if (((vlist[0] == v) && (vlist[1] == w)) || ((vlist[0] == w) && (vlist[1] == v)))
		id = vlist[2];
	else if (((vlist[1] == v) && (vlist[2] == w)) || ((vlist[1] == w) && (vlist[2] == v)))
		id = vlist[0];
	else
		id = vlist[1];
	return id;
}
void GeodesicPath::vertex_stepping(int v, int &w, int &f)
{
	w = m_mesh->vertex(v)->list()[0];
	double min = m_dist[w];
	for (int i = 1; i < m_mesh->vertex(v)->nNeighbor(); i++)
	{
		int x = m_mesh->vertex(v)->list()[i];
		if (m_dist[x] < min)
		{
			min = m_dist[x];
			w = x;
		}
	}
	int key1 = get_edge(v, w);
	int f1 = m_table[key1].face;
	int key2 = get_edge(w, v);
	int f2 = m_table[key2].face;
	int z1 = get_vertex_face(f1, v, w);
	int z2 = get_vertex_face(f2, v, w);
	f = (m_dist[z1] < m_dist[z2]) ? f1 : f2;
}

