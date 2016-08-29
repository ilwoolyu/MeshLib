/*************************************************
*	Slicer.h
*
*	Release: July 2011
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include "AABB.h"

#ifndef SLICER_HH_
#define SLICER_HH_

class Slicer: public AABB
{
private:
	class Point
	{
	public:
		MathVector pos;
		MathVector ID;
		MathVector nextID;
		bool trace;
		int faceID;
		Point *next;	// next / candidate (not yet connected)
		Point *prev;	// tail
		Point(void);
		Point(MathVector pos);
		bool operator <(const Point &p) const;
		bool operator >(const Point &p) const;
		bool operator ==(const Point &p) const;
		bool operator <(const MathVector &v) const;
		bool operator >(const MathVector &v) const;
		bool operator ==(const MathVector &v) const;
	};
	struct ContourHead
	{
		Point *begin;
		int size;
		bool closed;
	};
	std::vector<Point> m_slice;
	std::vector<ContourHead> m_group;

public:
	Slicer(void);
	Slicer(const Mesh *mesh);
	~Slicer(void);
	int slicing(const float a, const float b, const float c, const float d);
	void getSlice(float *plist, int group);
	const float *getEdgeID(int id, int group);
	int size(int group);

private:
	void searchTree(node *root, std::vector<int> *cand, const float a, const float b, const float c, const float d);
	bool boxIntersection(const float a, const float b, const float c, const float d, const float x0, const float x1, const float y0, const float y1, const float z0, const float z1);
	int contour(void);
	int tracePList(int headID, int id);
};

#endif
