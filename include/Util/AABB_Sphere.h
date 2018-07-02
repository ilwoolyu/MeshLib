/*************************************************
*	AABB_Sphere.h
*
*	Release: March 2017
*	Update: May 2018
*
*	Vanderbilt University
*	Electrical Engineering and Computer Science
*	
*	Ilwoo Lyu, ilwoo.lyu@vanderbilt.edu
*************************************************/

#include <vector>
#include "Mesh.h"
#include "Geom.h"

#ifndef AABB_SPHERE_HH_
#define AABB_SPHERE_HH_

using std::vector;

class AABB_Sphere
{
protected:
	struct node
	{
		node *left;
		node *right;
		float phi0, theta0;
		float phi1, theta1;
		vector<int> cand;
	};
	node *m_tree;
	const Mesh *m_mesh;
	bool *m_visited;
	float m_maxPhi;
	int m_pole0, m_pole1;
	float m_minTheta, m_maxTheta;

public:
	AABB_Sphere(void);
	AABB_Sphere(const Mesh *mesh);
	~AABB_Sphere(void);
	int closestFace(float *v, float *coeff);
	void update();
	void searchTree(const float *p, vector<int> *cand, const float eps1 = 0, const float eps2 = 0);

protected:
	void initTree(void);

private:
	int closestFace(float *v, float *coeff, vector<int> &xcand, const float eps1, const float eps2);
	int closestFaceBool(float *v, float *coeff, bool *xcand, const float eps1, const float eps2);
	void searchTree(const float phi, const float theta, node *root, vector<int> *cand, const float eps1 = 0, const float eps2 = 0);
	void searchTree(const float phi_min, const float theta_min, const float phi_max, const float theta_max, node *root, std::vector<int> *cand, float eps = 0, bool trace = false);
	node *construction(std::vector<float *> range, std::vector<int> cand);
	void boundingBox(const Face &f, float *r);
	void boundingBox(node *root);
	void updateTree(node *root);
	void deleteTree(node *root);
};

#endif
