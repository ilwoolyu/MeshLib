/*************************************************
*	AABB.h
*
*	Release: July 2011
*	Update: June 2020
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>
#include "Mesh.h"
#include "Geom.h"

#ifndef AABB_HH_
#define AABB_HH_

using std::vector;

class AABB
{
protected:
	struct node
	{
		node *left;
		node *right;
		float x0, y0, z0;
		float x1, y1, z1;
		std::vector<int> cand;
	};
	node *m_tree;
	const Mesh *m_mesh;
	bool *m_visited;

public:
	AABB(void);
	AABB(const Mesh *mesh);
	~AABB(void);
	int closestFace(float *v, float *coeff, float range = 0, float maxdist = 0);
	void update();
	void searchTree(const float *p, vector<int> *cand, float eps = 0);

protected:
	void initTree(void);

private:
	void searchTree(const float *p, node *root, vector<int> *cand, float eps = 0);
	void searchTree(float *pmin, float *pmax, node *root, std::vector<int> *cand, float eps = 0, bool trace = false);
	node *construction(std::vector<float *> range, std::vector<int> cand);
	void boundingBox(const Face &f, float *r);
	void boundingBox(node *root);
	void updateTree(node *root);
	void deleteTree(node *root);
};

#endif
