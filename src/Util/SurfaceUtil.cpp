/*************************************************
*	SurfaceUtil.cpp
*
*	Release: February 2015
*	Update: Jun 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <algorithm>
#include <cstring>
#include "Geom.h"
#include "SurfaceUtil.h"

void SurfaceUtil::smoothing(Mesh *mesh, int iter)
{
	for (int i = 0; i < iter; i++) smoothing(mesh);
}

void SurfaceUtil::smoothing(Mesh *mesh, float sigma, int iter)
{
	for (int i = 0; i < iter; i++) smoothing(mesh, sigma);
}

void SurfaceUtil::smoothing(const Mesh *mesh, int iter, float *data)
{
	for (int i = 0; i < iter; i++) smoothing(mesh, data);
}

void SurfaceUtil::smoothing(const Mesh *mesh, int iter, double *data)
{
	for (int i = 0; i < iter; i++) smoothing(mesh, data);
}

void SurfaceUtil::curvature(const Mesh *mesh, double *cmin, double *cmax, double **umin, double **umax, int nIterSmoothingTensor)
{
	const int n = mesh->nVertex();

	std::vector<double **> Tv;
	double *work = new double[9 * n];
	for (int i = 0; i < n; i++)
	{
		double **T = new double*[3];
		T[0] = &work[i * 9 + 0];
		T[1] = &work[i * 9 + 3];
		T[2] = &work[i * 9 + 6];
		Tv.push_back(T);
	}
	memset(work, 0, sizeof(double) * 9 * n);

	// edge-based tensor approximation
	tensor(mesh, Tv);

	// smoothing
	if (nIterSmoothingTensor > 0) smoothingTensor(mesh, Tv, nIterSmoothingTensor);
	
	double *work_min, *work_max;
	bool minNull = (umin == NULL), maxNull = (umax == NULL);
	if (minNull)
	{
		umin = new double*[n];
		work_min = new double[n * 3];
		for (int i = 0; i < n; i++)
			umin[i] = &work_min[i * 3];
	}
	if (maxNull)
	{
		umax = new double*[n];
		work_max = new double[n * 3];
		for (int i = 0; i < n; i++)
			umax[i] = &work_max[i * 3];
	}

	// principal curvature
	principalCurvature(Tv, n, cmin, cmax, umin, umax);

	for (int i = 0; i < n; i++)
		delete [] Tv[i];
	delete [] work;
	
	if (minNull)
	{
		delete [] work_min;
		delete [] umin;
		umin = NULL;
	}
	if (maxNull)
	{
		delete [] work_max;
		delete [] umax;
		umax = NULL;
	}
}

void SurfaceUtil::curvature(const Mesh *mesh, float *cmin, float *cmax, float **umin, float **umax, int nIterSmoothingTensor)
{
	const int n = mesh->nVertex();
	double *cmin_d = new double[n];
	double *cmax_d = new double[n];
	double **umin_d = new double*[n];
	double **umax_d = new double*[n];
	double *work_min = new double[n * 3];
	double *work_max = new double[n * 3];
	for (int i = 0; i < n; i++)
	{
		umin_d[i] = &work_min[i * 3];
		umax_d[i] = &work_max[i * 3];
	}
	curvature(mesh, cmin_d, cmax_d, umin_d, umax_d, nIterSmoothingTensor);
	for (int i = 0; i < n; i++)
	{
		cmin[i] = (float)cmin_d[i];
		cmax[i] = (float)cmax_d[i];
		for (int j = 0; j < 3; j++)
		{
			if (umin != NULL) umin[i][j] = (float)umin_d[i][j];
			if (umax != NULL) umax[i][j] = (float)umax_d[i][j];
		}
	}

	delete [] cmin_d;
	delete [] cmax_d;
	delete [] umin_d;
	delete [] umax_d;
	delete [] work_min;
	delete [] work_max;
}

void SurfaceUtil::sphere(Mesh *mesh, int type, int max_iter)
{
	const int n = mesh->nVertex();
	float **w = new float*[n];
	Vector *newv = new Vector[n];
	int c = 0;
	for (int i = 0; i < n; i++) c += mesh->vertex(i)->nNeighbor();
	float *work = new float[c]; memset(work, 0, sizeof(float) * c);
	c = 0;
	for (int i = 0; i < n; i++)
	{
		w[i] = &work[c];
		c += mesh->vertex(i)->nNeighbor();
	}
	mesh->centering();
	bool nInv;

	if (type == 0)
	{
		// area preservation
		for (int i = 0; i < mesh->nFace(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				int id = mesh->face(i)->list(j);
				int id1 = mesh->face(i)->list((j + 1) % 3);
				int id2 = mesh->face(i)->list((j + 2) % 3);
				Vector P = Vector(mesh->vertex(id)->fv(), mesh->vertex(id1)->fv());
				Vector Q = Vector(mesh->vertex(id)->fv(), mesh->vertex(id2)->fv());
				
				// angle
				/*float x = P.unit() * Q.unit();
				float weight = max(x / sqrt(1 - x * x), 1e-1);
				double weight = max(1 / tan(acos(x)), 1e-1);*/

				float weight = P.cross(Q).norm() / 2 / 12;
				for (int k = 0; k < mesh->vertex(id1)->nNeighbor(); k++)
					if (mesh->vertex(id1)->list(k) == id2) w[id1][k] += weight;
				for (int k = 0; k < mesh->vertex(id2)->nNeighbor(); k++)
					if (mesh->vertex(id2)->list(k) == id1) w[id2][k] += weight;
			}
		}
	}
	else
	{
		// angle preservation
		for (int i = 0; i < n; i++)
		{
			int nn = mesh->vertex(i)->nNeighbor();
			for (int j = 0; j < nn; j++)
			{
				int id1 = mesh->vertex(i)->list(j);
				int id2 = mesh->vertex(i)->list((j + 1) % nn);
				Vector P = Vector(mesh->vertex(i)->fv(), mesh->vertex(id1)->fv());
				Vector Q = Vector(mesh->vertex(i)->fv(), mesh->vertex(id2)->fv());

				double x = P.unit() * Q.unit();
				double weight = std::max(x / sqrt(1 - x * x), 1e-1);
				//double weight = max(1 / tan(acos(x)), 1e-1);

				for (int k = 0; k < mesh->vertex(id1)->nNeighbor(); k++)
					if (mesh->vertex(id1)->list(k) == id2) w[id1][k] += weight;
				for (int k = 0; k < mesh->vertex(id2)->nNeighbor(); k++)
					if (mesh->vertex(id2)->list(k) == id1) w[id2][k] += weight;
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		int nn = mesh->vertex(i)->nNeighbor();
		float sum = 0;
		for (int j = 0; j < nn; j++) sum += w[i][j];
		for (int j = 0; j < nn; j++) w[i][j] /= sum;
	}

	for (int i = 0; i < n; i++)
	{
		Vertex *v = (Vertex *)mesh->vertex(i);
		v->setVertex(Vector(v->fv()).unit().fv());
	}

	c = 0;
	do
	{
		/*mesh->updateNormal();
		for (int i = 0; i < n; i++)
			inv[i] = Vector(mesh->vertex(i)->fv()) * Vector(mesh->normal(i)->fv()) < 0;*/
		nInv = false;
		for (int i = 0; i < n; i++)
		{
			const int *neighbor = mesh->vertex(i)->list();
			const int nn = mesh->vertex(i)->nNeighbor();

			if (type == 0)
				newv[i] = mesh->vertex(i)->fv();	// area preservation
			else
				newv[i] = Vector({0, 0, 0});	// angle preservation
			
			for (int j = 0; j < nn; j++)
				newv[i] += Vector(mesh->vertex(neighbor[j])->fv()) * w[i][j];
			newv[i].unit();
		}
		for (int i = 0; i < n; i++)
		{
			Vertex *v = (Vertex *)mesh->vertex(i);
			v->setVertex(newv[i].fv());
		}
		for (int i = 0; i < mesh->nFace() && !nInv; i++)
		{
			Vector V0(mesh->vertex(mesh->face(i)->list(0))->fv());
			Vector V1(mesh->vertex(mesh->face(i)->list(1))->fv());
			Vector V2(mesh->vertex(mesh->face(i)->list(2))->fv());
			Vector V3 = V0 + V1 + V2;
			nInv = ((V1 - V0).cross(V2 - V0) * V3 < 0);
		}
		c++;
	} while (nInv && c < max_iter);

	c = 0;
	for (int i = 0; i < mesh->nFace(); i++)
	{
		Vector V0(mesh->vertex(mesh->face(i)->list(0))->fv());
		Vector V1(mesh->vertex(mesh->face(i)->list(1))->fv());
		Vector V2(mesh->vertex(mesh->face(i)->list(2))->fv());
		Vector V3 = V0 + V1 + V2;
		if ((V1 - V0).cross(V2 - V0) * V3 < 0) c++;
	}

	mesh->updateNormal();

	//delete [] inv;
	delete [] newv;
	delete [] work;
	delete [] w;
}

void SurfaceUtil::smoothing(Mesh *mesh)
{
	const int n = mesh->nVertex();
	Vector *newv = new Vector[n];

	for (int i = 0; i < n; i++)
	{
		const int *neighbor = mesh->vertex(i)->list();
		const int nn = mesh->vertex(i)->nNeighbor();
		
		for (int j = 0; j <  nn; j++)
		{
			newv[i] += Vector(mesh->vertex(neighbor[j])->fv());
		}
		newv[i] /= nn;
	}
	for (int i = 0; i < n; i++)
	{
		Vertex *v = (Vertex *)mesh->vertex(i);
		v->setVertex(newv[i].fv());
	}
	mesh->updateNormal();

	delete [] newv;
}

void SurfaceUtil::smoothing(Mesh *mesh, float sigma)
{
	const int n = mesh->nVertex();
	Vector *newv = new Vector[n];

	for (int i = 0; i < n; i++)
	{
		const int *neighbor = mesh->vertex(i)->list();
		const int nn = mesh->vertex(i)->nNeighbor();
		
		float weight = 0;
		for (int j = 0; j <  nn; j++)
		{
			float dist = Vector(mesh->vertex(neighbor[j])->fv(), mesh->vertex(i)->fv()).norm();
			float w = Statistics::normal_pdf(dist, 0, sigma);
			weight += w;
			newv[i] += Vector(mesh->vertex(neighbor[j])->fv()) * w;
		}
		newv[i] /= weight;
	}
	for (int i = 0; i < n; i++)
	{
		Vertex *v = (Vertex *)mesh->vertex(i);
		v->setVertex(newv[i].fv());
	}
	mesh->updateNormal();

	delete [] newv;
}

void SurfaceUtil::smoothing(const Mesh *mesh, float *data)
{
	const int n = mesh->nVertex();
	float *newv = new float[n];

	for (int i = 0; i < n; i++)
	{
		const int *neighbor = mesh->vertex(i)->list();
		const int nn = mesh->vertex(i)->nNeighbor();
		
		newv[i] = 0;
		for (int j = 0; j <  nn; j++)
		{
			newv[i] += data[neighbor[j]];
		}
		newv[i] /= nn;
	}

	memcpy(data, newv, sizeof(float) * n);

	delete [] newv;
}

void SurfaceUtil::smoothing(const Mesh *mesh, double *data)
{
	const int n = mesh->nVertex();
	double *newv = new double[n];

	for (int i = 0; i < n; i++)
	{
		const int *neighbor = mesh->vertex(i)->list();
		const int nn = mesh->vertex(i)->nNeighbor();
		
		newv[i] = 0;
		for (int j = 0; j <  nn; j++)
		{
			newv[i] += data[neighbor[j]];
		}
		newv[i] /= nn;
	}

	memcpy(data, newv, sizeof(double) * n);

	delete [] newv;
}

void SurfaceUtil::tensor(const Mesh *mesh, std::vector<double **> &Tv)
{
	const int n = mesh->nVertex();
	bool *check = new bool[n];	// flag per each vertex
	memset(check, 0, sizeof(bool) * n);

	// get an edge list
	std::vector<edge> edgeList;
	for (int i = 0; i < n; i++)
	{
		// check its neighbroing points
		for (int j = 0; j < mesh->vertex(i)->nNeighbor(); j++)
		{
			int id = mesh->vertex(i)->list(j);

			if (!check[id])
			{
				edge e;
				e.vid1 = i;
				e.vid2 = id;
				e.fid1 = -1;
				e.fid2 = -1;
				e.tan = Vector(mesh->vertex(e.vid1)->fv(), mesh->vertex(e.vid2)->fv());
				e.len = e.tan.norm();
				e.tan.unit();
				if (i > id) std::swap(e.vid1, e.vid2);
				edgeList.push_back(e);
			}
		}
		check[i] = true;
	}

	std::sort(edgeList.begin(), edgeList.end()); 

	for (int i = 0; i < mesh->nFace(); i++)
	{
		int id[3] = {mesh->face(i)->list(0), mesh->face(i)->list(1), mesh->face(i)->list(2)};
		edge e;
		std::vector<edge>::iterator it;

		// vertex 0 and 1
		e.vid1 = id[0]; e.vid2 = id[1];
		if (id[0] > id[1]) std::swap(e.vid1, e.vid2);
		it = std::lower_bound(edgeList.begin(), edgeList.end(), e);
		if (id[0] < id[1]) (*it).fid1 = i;
		else (*it).fid2 = i;

		// vertex 1 and 2
		e.vid1 = id[1]; e.vid2 = id[2];
		if (id[1] > id[2]) std::swap(e.vid1, e.vid2);
		it = std::lower_bound(edgeList.begin(), edgeList.end(), e);
		if (id[1] < id[2]) (*it).fid1 = i;
		else (*it).fid2 = i;

		// vertex 0 and 1
		e.vid1 = id[2]; e.vid2 = id[0];
		if (id[2] > id[0]) std::swap(e.vid1, e.vid2);
		it = std::lower_bound(edgeList.begin(), edgeList.end(), e);
		if (id[2] < id[0]) (*it).fid1 = i;
		else (*it).fid2 = i;
	}

	double meanLen = 0;
	for (int i = 0; i < edgeList.size(); i++)
		meanLen += edgeList[i].len;
	meanLen /= edgeList.size();

	// tensor
	for (int i = 0; i < edgeList.size(); i++)
	{
		int fid1 = edgeList[i].fid1;
		int fid2 = edgeList[i].fid2;
		if (fid1 == -1 || fid2 == -1) continue;
		
		Vector fn1(mesh->face(fid1)->faceNormal().fv()), fn2(mesh->face(fid2)->faceNormal().fv());
		double beta = fn1 * fn2;
		Vector cp = fn1.cross(fn2);

		if (beta > 1) beta = 1;
		else if (beta < -1) beta = -1;
		beta = (cp * edgeList[i].tan > 0) ? acos(beta): -acos(beta);	// convex/concave

		// tensor
		double bLen = beta * edgeList[i].len / meanLen;
		for (int j = 0; j < 3; j++)
			for (int k = j; k < 3; k++)
				edgeList[i].T[j][k] = edgeList[i].T[k][j] = edgeList[i].tan[j] * edgeList[i].tan[k] * bLen;

		// propagation over the two vertices
		for (int j = 0; j < 3; j++)
		{
			for (int k = j; k < 3; k++)
			{
				Tv[edgeList[i].vid1][j][k] += edgeList[i].T[j][k] / (double)mesh->vertex(edgeList[i].vid1)->nNeighbor();
				Tv[edgeList[i].vid1][k][j] = Tv[edgeList[i].vid1][j][k];
				
				Tv[edgeList[i].vid2][j][k] += edgeList[i].T[j][k] / (double)mesh->vertex(edgeList[i].vid2)->nNeighbor();
				Tv[edgeList[i].vid2][k][j] = Tv[edgeList[i].vid2][j][k];
			}
		}
	}
	delete [] check;
}

void SurfaceUtil::smoothingTensor(const Mesh *mesh, const std::vector<double **> Tv, int iter)
{
	const int n = mesh->nVertex();
	std::vector<double **> TvT;
	double *workT = new double[9 * n];
	for (int i = 0; i < n; i++)
	{
		double **T = new double*[3];
		T[0] = &workT[i * 9 + 0];
		T[1] = &workT[i * 9 + 3];
		T[2] = &workT[i * 9 + 6];
		TvT.push_back(T);
	}
	for (int s = 0; s < iter; s++)
	{
		memset(workT, 0, sizeof(double) * 9 * n);
		for (int i = 0; i < n; i++)
		{
			double weight = 1;
			for (int j = 0; j < 3; j++)
			{
				for (int k = j; k < 3; k++)
				{
					TvT[i][j][k] += Tv[i][j][k];
					TvT[i][k][j] = TvT[i][j][k];
				}
			}
			//double weight = 0;
			for (int nn = 0; nn < mesh->vertex(i)->nNeighbor(); nn++)
			{
				int id = mesh->vertex(i)->list(nn);
				//double w = 1.0f / Vector(mesh->vertex(i)->fv(), mesh->vertex(id)->fv()).norm();
				double w = 1.0f;
				for (int j = 0; j < 3; j++)
				{
					for (int k = j; k < 3; k++)
					{
						TvT[i][j][k] += Tv[id][j][k] * w;
						TvT[i][k][j] = TvT[i][j][k];
					}
				}
				weight += w;
			}
			for (int j = 0; j < 3; j++)
			{
				for (int k = j; k < 3; k++)
				{
					TvT[i][j][k] /= weight;
					TvT[i][k][j] = TvT[i][j][k];
				}
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = j; k < 3; k++)
				{
					Tv[i][j][k] = TvT[i][j][k];
					Tv[i][k][j] = Tv[i][j][k];
				}
			}
		}
	}
	
	for (int i = 0; i < n; i++)
		delete [] TvT[i];
	delete [] workT;
}

void SurfaceUtil::principalCurvature(const std::vector<double **> Tv, const int n, double *cmin, double *cmax, double **umin, double **umax)
{
	for (int i = 0; i < n; i++)
	{
		double lambda[3];
		double eigv[3][3];
		double T[3][3] = {{Tv[i][0][0], Tv[i][0][1], Tv[i][0][2]}, {Tv[i][1][0], Tv[i][1][1], Tv[i][1][2]}, {Tv[i][2][0], Tv[i][2][1], Tv[i][2][2]}};

		LinearAlgebra::eig3symmetric(T, lambda, eigv);

		if ((fabs(lambda[0]) > fabs(lambda[1]) && fabs(lambda[0]) < fabs(lambda[2])) || 
			(fabs(lambda[0]) < fabs(lambda[1]) && fabs(lambda[0]) > fabs(lambda[2])))
		{
			cmin[i] = lambda[0];
			memcpy(umax[i], eigv[0], sizeof(double) * 3);
		}
		else if ((fabs(lambda[1]) > fabs(lambda[2]) && fabs(lambda[1]) < fabs(lambda[0])) ||
				(fabs(lambda[1]) < fabs(lambda[2]) && fabs(lambda[1]) > fabs(lambda[0])))
		{
			cmin[i] = lambda[1];
			memcpy(umax[i], eigv[1], sizeof(double) * 3);
		}
		else if ((fabs(lambda[2]) > fabs(lambda[0]) && fabs(lambda[2]) < fabs(lambda[1])) ||
				(fabs(lambda[2]) < fabs(lambda[0]) && fabs(lambda[2]) > fabs(lambda[1])))
		{
			cmin[i] = lambda[2];
			memcpy(umax[i], eigv[2], sizeof(double) * 3);
		}

		if (fabs(lambda[0]) > fabs(lambda[1]) && fabs(lambda[0]) > fabs(lambda[2]))
		{
			cmax[i] = lambda[0];
			memcpy(umin[i], eigv[0], sizeof(double) * 3);
		}
		else if (fabs(lambda[1]) > fabs(lambda[2]) && fabs(lambda[1]) > fabs(lambda[0]))
		{
			cmax[i] = lambda[1];
			memcpy(umin[i], eigv[1], sizeof(double) * 3);
		}
		else if (fabs(lambda[2]) > fabs(lambda[0]) && fabs(lambda[2]) > fabs(lambda[1]))
		{
			cmax[i] = lambda[2];
			memcpy(umin[i], eigv[2], sizeof(double) * 3);
		}

		// swap
		if (cmax[i] < cmin[i])
		{
			std::swap(cmax[i], cmin[i]);
			std::swap(umax[i], umin[i]);
		}
	}
}

