/*************************************************
*	Slicer.cpp
*
*	Release: July 2011
*	Update: April 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <algorithm>
#include <cstring>
#include "Slicer.h"

using namespace std;

Slicer::Point::Point(void)
{
	this->pos = 0.0f;
	next = NULL;
	prev = NULL;
	trace = false;
}
Slicer::Point::Point(MathVector pos)
{
	this->pos = pos;
	next = NULL;
	prev = NULL;
	trace = false;
}
bool Slicer::Point::operator <(const Point &p) const
{
	return this->ID < p.ID;
}
bool Slicer::Point::operator >(const Point &p) const
{
	return this->ID > p.ID;
}
bool Slicer::Point::operator ==(const Point &p) const
{
	return this->ID == p.ID;
}
bool Slicer::Point::operator <(const MathVector &v) const
{
	return this->ID < v;
}
bool Slicer::Point::operator >(const MathVector &v) const
{
	return this->ID > v;
}
bool Slicer::Point::operator ==(const MathVector &v) const
{
	return this->ID == v;
}

Slicer::Slicer(void)
{
	m_mesh = NULL;
	m_slice.clear();
}
Slicer::Slicer(const Mesh *mesh)
{
	m_mesh = mesh;
//	printf("AABB Tree Construction.. ");
	initTree();
	update();
//	printf("Done\n");
}
Slicer::~Slicer(void)
{
	m_slice.clear();
}
int Slicer::slicing(const float a, const float b, const float c, const float d)
{
	vector<int> cand;
	searchTree(m_tree, &cand, a, b, c, d);

	sort(cand.begin(), cand.end());
	cand.erase(unique(cand.begin(), cand.end()), cand.end());

	m_slice.clear();
	m_group.clear();

	for (int i = 0; i < cand.size(); i++)
	{
		const float *t0 = m_mesh->face(cand[i])->vertex(0)->fv();
		const float *t1 = m_mesh->face(cand[i])->vertex(1)->fv();
		const float *t2 = m_mesh->face(cand[i])->vertex(2)->fv();

		float p0[3], p1[3];
		int nv = Coordinate::intersection(t0, t1, t2, a, b, c, d, p0, p1);
		if (nv == 2)
		{
			// consistent direction
			MathVector N = m_mesh->face(cand[i])->faceNormal().fv();
			MathVector P0 = p0;
			MathVector P1 = p1;
			MathVector P = P0 + N;
			float p[3];
			Coordinate::proj2plane(a, b, c, d, P(), p);
			P = p;

			// edge id
			float it0 = a * t0[0] + b * t0[1] + c * t0[2] + d;
			float it1 = a * t1[0] + b * t1[1] + c * t1[2] + d;
			float it2 = a * t2[0] + b * t2[1] + c * t2[2] + d;

			/*float eps = 1e-8;
			if (fabs(it0) < eps) it0 = 0;
			if (fabs(it1) < eps) it1 = 0;
			if (fabs(it2) < eps) it2 = 0;*/

			if (MathVector(t0) == MathVector(p0) || MathVector(t0) == MathVector(p1))
				it0 = 0;
			if (MathVector(t1) == MathVector(p0) || MathVector(t1) == MathVector(p1))
				it1 = 0;
			if (MathVector(t2) == MathVector(p0) || MathVector(t2) == MathVector(p1))
				it2 = 0;

			int id0[3] = {0, 0, 0}, id1[3] = {0, 0, 0};
			//if (fabs(it1) < eps && fabs(it2) < eps)
			if (fabs(it1) == 0 && fabs(it2) == 0)
			{
				id0[0] = id0[1] = m_mesh->face(cand[i])->vertex(1)->id();
				id1[0] = id1[1] = m_mesh->face(cand[i])->vertex(2)->id();
			}
			//else if (fabs(it2) < eps && fabs(it0) < eps)
			else if (fabs(it2) == 0 && fabs(it0) == 0)
			{
				id0[0] = id0[1] = m_mesh->face(cand[i])->vertex(2)->id();
				id1[0] = id1[1] = m_mesh->face(cand[i])->vertex(0)->id();
			}
			//else if (fabs(it0) < eps && fabs(it1) < eps)
			else if (fabs(it0) == 0 && fabs(it1) == 0)
			{
				id0[0] = id0[1] = m_mesh->face(cand[i])->vertex(0)->id();
				id1[0] = id1[1] = m_mesh->face(cand[i])->vertex(1)->id();
			}
			else if (it1 * it2 >= 0)
			{
				//if (fabs(it1) < eps)
				if (fabs(it1) == 0)
					id0[0] = m_mesh->face(cand[i])->vertex(1)->id();
				else
					id0[0] = m_mesh->face(cand[i])->vertex(0)->id();
				//if (fabs(it2) < eps)
				if (fabs(it2) == 0)
					id1[0] = m_mesh->face(cand[i])->vertex(2)->id();
				else
					id1[0] = m_mesh->face(cand[i])->vertex(0)->id();
				id0[1] = m_mesh->face(cand[i])->vertex(1)->id();
				id1[1] = m_mesh->face(cand[i])->vertex(2)->id();
			}
			else if (it2 * it0 >= 0)
			{
				//if (fabs(it2) < eps)
				if (fabs(it2) == 0)
					id0[0] = m_mesh->face(cand[i])->vertex(2)->id();
				else
					id0[0] = m_mesh->face(cand[i])->vertex(1)->id();
				//if (fabs(it0) < eps)
				if (fabs(it0) == 0)
					id1[0] = m_mesh->face(cand[i])->vertex(0)->id();
				else
					id1[0] = m_mesh->face(cand[i])->vertex(1)->id();
				id0[1] = m_mesh->face(cand[i])->vertex(2)->id();
				id1[1] = m_mesh->face(cand[i])->vertex(0)->id();
			}
			else if (it0 * it1 >= 0)
			{
				//if (fabs(it0) < eps)
				if (fabs(it0) == 0)
					id0[0] = m_mesh->face(cand[i])->vertex(0)->id();
				else
					id0[0] = m_mesh->face(cand[i])->vertex(2)->id();
				//if (fabs(it1) < eps)
				if (fabs(it1) == 0)
					id1[0] = m_mesh->face(cand[i])->vertex(1)->id();
				else
					id1[0] = m_mesh->face(cand[i])->vertex(2)->id();
				id0[1] = m_mesh->face(cand[i])->vertex(0)->id();
				id1[1] = m_mesh->face(cand[i])->vertex(1)->id();
			}
			if (id0[0] > id0[1]) swap(id0[0], id0[1]);
			if (id1[0] > id1[1]) swap(id1[0], id1[1]);

			MathVector PP0 = P0 - P, PP1 = P1 - P;
			MathVector PCP = PP0.cross(PP1);

			if (PCP * MathVector(a, b, c) > 0)
			{
				swap(p0, p1);
				swap(id0, id1);
			}

			// new item
			Point item;
			item.pos = p0;
			item.ID = id0;
			item.faceID = cand[i];
			item.nextID = id1;
			m_slice.push_back(item);
		}
	}
	sort(m_slice.begin(), m_slice.end());

	return contour();
}
void Slicer::getSlice(float *plist, int group)
{
	Point *iter = m_group[group].begin;
	for (int i = 0; i < m_group[group].size; i++)
	{
		memcpy(&plist[i * 3], iter->pos.fv(), sizeof(float) * 3);
		iter = iter->next;
	}
}
void Slicer::getSliceFaceIdx(int *plist, int group)
{
	Point *iter = m_group[group].begin;
	for (int i = 0; i < m_group[group].size; i++)
	{
		plist[i] = iter->faceID;
		iter = iter->next;
	}
}
const float * Slicer::getEdgeID(int id, int group)
{
	Point *iter = m_group[group].begin;
	for (int i = 0; i < m_group[group].size; i++)
	{
		if (i == id) return iter->ID.fv();
		iter = iter->next;
	}
	return NULL;
}
int Slicer::size(int group)
{
	return m_group[group].size;
}
void Slicer::searchTree(node *root, vector<int> *cand, const float a, const float b, const float c, const float d)
{
	if (root->left != NULL && 
		boxIntersection(a, b, c, d, root->left->x0, root->left->x1, root->left->y0, root->left->y1, root->left->z0, root->left->z1))
		searchTree(root->left, cand, a, b, c, d);

	if (root->right != NULL && 
		boxIntersection(a, b, c, d, root->right->x0, root->right->x1, root->right->y0, root->right->y1, root->right->z0, root->right->z1))
		searchTree(root->right, cand, a, b, c, d);

	if (root->left == NULL && root->right == NULL)
		for (int i = 0; i < root->cand.size(); i++)
			cand->push_back(root->cand[i]);
}
bool Slicer::boxIntersection(const float a, const float b, const float c, const float d, const float x0, const float x1, const float y0, const float y1, const float z0, const float z1)
{
	float i1 = (a * x0 + b * y0 + c * z0 + d);
	float i2 = (a * x1 + b * y0 + c * z0 + d);
	float i3 = (a * x0 + b * y1 + c * z0 + d);
	float i4 = (a * x1 + b * y1 + c * z0 + d);
	float i5 = (a * x0 + b * y0 + c * z1 + d);
	float i6 = (a * x1 + b * y0 + c * z1 + d);
	float i7 = (a * x0 + b * y1 + c * z1 + d);
	float i8 = (a * x1 + b * y1 + c * z1 + d);

	if ((i1 > 0 && i2 > 0 && i3 > 0 && i4 > 0 && i5 > 0 && i6 > 0 && i7 > 0 && i8 > 0) ||
		(i1 < 0 && i2 < 0 && i3 < 0 && i4 < 0 && i5 < 0 && i6 < 0 && i7 < 0 && i8 < 0))
		return false;
	else
		return true;
}
int Slicer::contour(void)
{
	bool closed = false;
	for (int i = 0; i < m_slice.size(); i++)
	//for (int i = m_slice.size() - 1; i >= 0; i--)
	{
		if (!m_slice[i].trace)
		{
			int n = tracePList(i, i);
			if (n > 0)
			{
				ContourHead ch;
				ch.begin = &m_slice[i];
				ch.size = n;
				ch.closed = (ch.begin == ch.begin->prev->next);

				if (!ch.closed)
				{
					// partial loop
					bool inside = false;
					Point *part = ch.begin->prev->next;
					Point *iter = ch.begin;

					// check if the next candidate is inside the current loop
					for (int j = 0; j < ch.size; j++, iter = iter->next)
					{
						if (iter == part)
						{
							inside = true;
							break;
						}
					}
					if (inside)
					{
						ContourHead chLoop;
						swap(ch.begin->prev, part->prev);
						
						// previous group
						for (n = 0, iter = ch.begin; iter != part; iter = iter->next, n++);
						ch.size = n;

						// new group
						chLoop.begin = part;
						chLoop.closed = true;
						for (n = 1, iter = chLoop.begin; iter != chLoop.begin->prev; iter = iter->next, n++);
						chLoop.size = n;
						m_group.push_back(chLoop);
					}
				}

				closed = closed || !ch.closed;
				m_group.push_back(ch);
			}
		}
	}

	// connected components
	int gSize = m_group.size();
	while (closed)
	{
		int i = 0;
		while (i < m_group.size())
		{
			//cout << "Group " << i << ": " << m_group[i].closed << endl;
			if (!m_group[i].closed)
			{
				for (int c = 0; c < m_group.size(); c++)
				{
					if (m_group[c].closed) continue;
					Point *tail = m_group[i].begin->prev->next;
					Point *head = m_group[c].begin;
					if (tail == head)
					{
						swap(m_group[i].begin->prev, m_group[c].begin->prev);
						m_group[i].size += m_group[c].size;
						m_group.erase(m_group.begin() + c);
						break;
					}
				}
			}
			i++;
		}
		closed = (gSize != m_group.size());
		gSize = m_group.size();
	}

	return m_group.size();
}
int Slicer::tracePList(int headID, int id)
{
	// next target
	int nextID = lower_bound(m_slice.begin(), m_slice.end(), m_slice[id].nextID) - m_slice.begin();
	/*int nextID;
	for (int i = 0; i < m_slice.size(); i++)
	{
		if (m_slice[i].ID == m_slice[id].nextID && i != id)
		{
			nextID = i;
			break;
		}
	}*/

	m_slice[id].trace = true;
	m_slice[id].next = &m_slice[nextID];

	if (!m_slice[nextID].trace)
	{
		m_slice[nextID].prev = &m_slice[id];

		return tracePList(headID, nextID) + 1;
	}
	/*else if (m_slice[nextID].faceID != headFaceID)
	{
		cout << m_slice[nextID].trace << endl;

		int newNextID = m_slice.size();
		Point *iter;
		int n = 0;

		for (iter = &m_slice[nextID]; iter->trace; iter = iter->prev, n++)
		{
			Point item = m_slice[nextID];
			item.trace = false;
			item.next = &item;
			m_slice.push_back(item);
		}
		m_slice[nextID].prev = &m_slice[id];
		m_slice[id].next = &m_slice[newNextID];

		sort(m_slice.begin(), m_slice.end());
		
		int nextID = lower_bound(m_slice.begin(), m_slice.end(), iter->nextID) - m_slice.begin();

		return tracePList(headFaceID, nextID) + n;
	}*/
	else
	{
		m_slice[headID].prev = &m_slice[id];

		int nextID = -1;
		for (int i = 0; i < m_slice.size(); i++)
		{
			if (m_slice[i].ID == m_slice[id].nextID && i != id)
			{
				nextID = i;
				break;
			}
		}

		return 1;
	}
}

