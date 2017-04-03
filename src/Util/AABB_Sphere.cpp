/*************************************************
*	AABB_Sphere.cpp
*
*	Release: March 2017
*	Update: April 2017
*
*	Vanderbilt University
*	Electrical Engineering and Computer Science
*	
*	Ilwoo Lyu, ilwoo.lyu@vanderbilt.edu
*************************************************/

#include <stdio.h>
#include <algorithm>
#include <float.h>
#include <cstring>
#include "AABB_Sphere.h"

using namespace std;

AABB_Sphere::AABB_Sphere(void)
{
}
AABB_Sphere::AABB_Sphere(const Mesh *mesh)
{
	m_mesh = mesh;
	initTree();
	update();
}
AABB_Sphere::~AABB_Sphere(void)
{
}
int AABB_Sphere::closestFace(float *v, float *coeff)
{
	vector<int> xcand;

	int index = -1;
	float rad = 0.0174533; // 0.0174533 rad = 1 degree
	for (int trial = 0; index == -1 && trial <= 9; trial++)
	{
		index = closestFace(v, coeff, xcand, rad * trial);
		rad *= 2;
	}
	return index;
}
int AABB_Sphere::closestFace(float *v, float *coeff, vector<int> &xcand, const float eps)
{
	int index = -1;

	float phi, theta;
	Coordinate::cart2sph(v, &phi, &theta);

	vector<int> cand_all;
	searchTree(phi, theta, m_tree, &cand_all, eps);
	if (phi + 2 * PI <= m_maxPhi) searchTree(phi + 2 * PI, theta, m_tree, &cand_all, eps);
	
	if (cand_all.empty()) return -1;

	sort(cand_all.begin(), cand_all.end());
	cand_all.erase(unique(cand_all.begin(), cand_all.end()), cand_all.end());

	vector<int> cand;	// only consider unseen triangles
	set_difference(cand_all.begin(), cand_all.end(), xcand.begin(), xcand.end(), inserter(cand, cand.begin()));

	float min_dist = FLT_MAX;
	for (int i = 0; i < cand.size(); i++)
	{
		const float *a = m_mesh->face(cand[i])->vertex(0)->fv();
		const float *b = m_mesh->face(cand[i])->vertex(1)->fv();
		const float *c = m_mesh->face(cand[i])->vertex(2)->fv();

		float tcoeff[3];
		// projection of v onto the triangle
		float v_proj[3];
		Vector N = Vector(a, b).cross(Vector(b, c)).unit();
		const float *n = N.fv();
		Vector V_proj = Vector(v) * (Vector(a) * N) / (Vector(v) * N);	// scaling a query vector
		
		Coordinate::cart2bary((float *)a, (float *)b, (float *)c, (float *)V_proj.fv(), tcoeff);
		if (tcoeff[0] >= 0 && tcoeff[1] >= 0 && tcoeff[2] >= 0)
		{
			// closest distance
			float dist = Coordinate::dpoint2tri(a, b, c, v);
			if (dist < min_dist)
			{
				index = cand[i];
				min_dist = dist;
			
				// bary centric
				coeff[0] = tcoeff[0];
				coeff[1] = tcoeff[1];
				coeff[2] = tcoeff[2];
			}
		}
	}
	xcand.insert(xcand.end(), cand.begin(), cand.end());
	
	return index;
}
void AABB_Sphere::update()
{
	updateTree(m_tree);
}
void AABB_Sphere::initTree(void)
{
	vector<int> cand;
	vector<float *> range;
	for (int i = 0; i < m_mesh->nFace(); i++)
	{
		const Face &f = *m_mesh->face(i);
		float *r = new float[4];
		boundingBox(f, r);
		range.push_back(r);
	}
	for (int i = 0; i < m_mesh->nFace(); i++) cand.push_back(i);

	m_tree = construction(range, cand);
}
void AABB_Sphere::searchTree(const float *p, vector<int> *cand, float eps)
{
	float phi, theta;
	Coordinate::cart2sph(p, &phi, &theta);

	searchTree(phi, theta, m_tree, cand, eps);
	if (phi + 2 * PI <= m_maxPhi) searchTree(phi + 2 * PI, theta, m_tree, cand, eps);	
}
void AABB_Sphere::searchTree(const float phi, const float theta, node *root, vector<int> *cand, float eps)
{
	if (root->left != NULL && 
		phi >= root->left->phi0 - eps && theta >= root->left->theta0 - eps &&
		phi <= root->left->phi1 + eps && theta <= root->left->theta1 + eps)
		searchTree(phi, theta, root->left, cand, eps);

	if (root->right != NULL &&
		phi >= root->right->phi0 - eps && theta >= root->right->theta0 - eps &&
		phi <= root->right->phi1 + eps && theta <= root->right->theta1 + eps)
		searchTree(phi, theta, root->right, cand, eps);

	if (root->left == NULL && root->right == NULL)
		for (int i = 0; i < root->cand.size(); i++)
			cand->push_back(root->cand[i]);
}
AABB_Sphere::node *AABB_Sphere::construction(const vector<float *> range, const vector<int> cand)
{
	vector<float *> lrange, rrange;
	vector<int> left, right;
	vector<float *> lrrange;
	vector<int> lr;

	if (cand.empty())
	{
		return NULL;
	}
	else if (cand.size() == 1)
	{
		node *elem = new node();
		elem->left = NULL;
		elem->right = NULL;
		elem->cand = cand;

		return elem;
	}
	else if (cand.size() == 2)
	{
		node *elem = new node();
		lrange.push_back(range[0]);
		left.push_back(cand[0]);
		rrange.push_back(range[1]);
		right.push_back(cand[1]);
		elem->left = construction(lrange, left);
		elem->right = construction(rrange, right);
		
		return elem;
	}
	else
	{
		// maximum length
		float r[4];
		for (int i = 0; i < 4; i++) r[i] = range[0][i];

		for (int i = 1; i < range.size(); i++)
		{
			for (int dim = 0; dim < 2; dim++)
			{
				if (r[dim * 2 + 0] > range[i][dim * 2 + 0]) r[dim * 2 + 0] = range[i][dim * 2 + 0];
				if (r[dim * 2 + 1] < range[i][dim * 2 + 1]) r[dim * 2 + 1] = range[i][dim * 2 + 1];
			}
		}

		// find the best divider
		int cdim = 0;
		float cthreshold = 0;
		float score = -1;
		for (int dim = 0; dim < 2; dim++)
		{
			// median
			std::vector<float> data;
			for (int i = 0; i < range.size(); i++)
			{
				data.push_back(range[i][dim * 2 + 0]);
				data.push_back(range[i][dim * 2 + 1]);
			}

			// sort for median
			std::sort(data.begin(), data.end());
			float threshold = data[(int)(data.size() / 2)];

			int na = 0, nb = 0, nc = 0;
			for (int i = 0; i < range.size(); i++)
			{
				float min = range[i][dim * 2 + 0];
				float max = range[i][dim * 2 + 1];
				if (threshold >= max) na++;
				else if (threshold <= min) nb++;
				else nc++;
			}

			if ((na > 0 && nb > 0) || nc > 0)
			{
				if (score < (float)na * na + (float)nb * nb)
				{
					cdim = dim;
					cthreshold = threshold;
					score = na * na + nb * nb;
				}
			}
		}

		// classification
		for (int i = 0; i < range.size(); i++)
		{
			float min = range[i][cdim * 2 + 0];
			float max = range[i][cdim * 2 + 1];
			if (cthreshold >= max)
			{
				left.push_back(cand[i]);
				lrange.push_back(range[i]);
			}
			else if (cthreshold <= min)
			{
				right.push_back(cand[i]);
				rrange.push_back(range[i]);
			}
			else
			{
				lr.push_back(cand[i]);
				lrrange.push_back(range[i]);
			}
		}
	}
	// new node
	node *elem = NULL;
	node *lnode, *rnode;

	// even distribution
	if (!lr.empty())
	{
		for (int i = 0; i < lr.size(); i++)
		{
			if (left.size() < right.size())
			{
				left.push_back(lr[i]);
				lrange.push_back(lrrange[i]);
			}
			else
			{
				right.push_back(lr[i]);
				rrange.push_back(lrrange[i]);
			}
		}
	}

	if (left.empty() || right.empty())
	{
		elem = new node();
		elem->left = NULL;
		elem->right = NULL;
		elem->cand = cand;
	}
	else
	{
		lnode = construction(lrange, left);
		rnode = construction(rrange, right);

		if (lnode != NULL && rnode == NULL)
		{
			elem = lnode;
		}
		else if (lnode == NULL && rnode != NULL)
		{
			elem = rnode;
		}
		else
		{
			elem = new node();
			elem->left = lnode;
			elem->right = rnode;
			if (lnode == NULL && rnode == NULL)
				elem->cand = cand;
		}
	}
	return elem;
}
void AABB_Sphere::boundingBox(const Face &f, float *r)
{
	float minphi = 2 * PI, maxphi = -2 * PI;
	float mintheta = 2 * PI, maxtheta = -2 * PI;

	for (int i = 0; i < 3; i++)
	{
		const Vertex *v = f.vertex(i);
		float phi, theta;
		Coordinate::cart2sph(v->fv(), &phi, &theta);

		if (minphi > phi) minphi = phi;
		if (maxphi < phi) maxphi = phi;
		if (mintheta > theta) mintheta = theta;
		if (maxtheta < theta) maxtheta = theta;
	}
	if (maxphi - minphi > PI)
	{
		minphi += 2 * PI;
		swap(minphi, maxphi);
		if (m_maxPhi < maxphi) m_maxPhi = maxphi;
	}
	r[0] = minphi; r[1] = maxphi;
	r[2] = mintheta; r[3] = maxtheta;
}
void AABB_Sphere::boundingBox(node *root)
{
	if (root->left == NULL && root->right == NULL)
	{
		const Face &f = *m_mesh->face(root->cand[0]);
		float r[4];
		boundingBox(f, r);
		root->phi0 = r[0]; root->theta0 = r[2];
		root->phi1 = r[1]; root->theta1 = r[3];
		for (int i = 1; i < root->cand.size(); i++)
		{
			Face f = *m_mesh->face(root->cand[i]);
			float r[4];
			boundingBox(f, r);
			root->phi0 = (root->phi0 < r[0]) ? root->phi0: r[0];
			root->theta0 = (root->theta0 < r[2]) ? root->theta0: r[2];
			root->phi1 = (root->phi1 > r[1]) ? root->phi1: r[1];
			root->theta1 = (root->theta1 > r[3]) ? root->theta1: r[3];
		}
	}
	else if (root->left != NULL && root->right == NULL)
	{
		root->phi0 = root->left->phi0;
		root->theta0 = root->left->theta0;
		root->phi1 = root->left->phi1;
		root->theta1 = root->left->theta1;
	}
	else if (root->left == NULL && root->right != NULL)
	{
		root->phi0 = root->right->phi0;
		root->theta0 = root->right->theta0;
		root->phi1 = root->right->phi1;
		root->theta1 = root->right->theta1;
	}
	else
	{
		//assert(root->left != NULL && root->right != NULL);
		root->phi0 = (root->left->phi0 < root->right->phi0) ? root->left->phi0: root->right->phi0;
		root->theta0 = (root->left->theta0 < root->right->theta0) ? root->left->theta0: root->right->theta0;
		root->phi1 = (root->left->phi1 > root->right->phi1) ? root->left->phi1: root->right->phi1;
		root->theta1 = (root->left->theta1 > root->right->theta1) ? root->left->theta1: root->right->theta1;
	}
}
void AABB_Sphere::updateTree(node *root)
{
	if (root->left != NULL) updateTree(root->left);
	if (root->right != NULL) updateTree(root->right);
	
	m_maxPhi = 2 * PI;
	boundingBox(root);
}

