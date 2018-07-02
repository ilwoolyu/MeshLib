/*************************************************
*	Geom.cpp
*
*	Release: July 2011
*	Update: October 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <cstring>
#include "Geom.h"

MathVector::MathVector(void)
{
	m_vector[0] = m_vector[1] = m_vector[2] = 0;
}

MathVector::MathVector(const float *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
}

MathVector::MathVector(const float *v1, const float *v2)
{
	m_vector[0] = v2[0] - v1[0];
	m_vector[1] = v2[1] - v1[1];
	m_vector[2] = v2[2] - v1[2];
}

MathVector::MathVector(const float v1, const float v2, const float v3)
{
	m_vector[0] = v1;
	m_vector[1] = v2;
	m_vector[2] = v3;
}

MathVector::MathVector(const int v1, const int v2, const int v3)
{
	m_vector[0] = (const float)v1;
	m_vector[1] = (const float)v2;
	m_vector[2] = (const float)v3;
}

MathVector::MathVector(const MathVector &v)
{
	m_vector[0] = v.m_vector[0];
	m_vector[1] = v.m_vector[1];
	m_vector[2] = v.m_vector[2];
}

MathVector::~MathVector(void)
{
}

MathVector MathVector::cross(const MathVector v)
{
	MathVector p;
	p.m_vector[0] = m_vector[1] * v.m_vector[2] - m_vector[2] * v.m_vector[1];
	p.m_vector[1] = m_vector[2] * v.m_vector[0] - m_vector[0] * v.m_vector[2];
	p.m_vector[2] = m_vector[0] * v.m_vector[1] - m_vector[1] * v.m_vector[0];
	return p;
}

float MathVector::norm(void)
{
	return sqrt(m_vector[0] * m_vector[0] + m_vector[1] * m_vector[1] + m_vector[2] * m_vector[2]);
}

MathVector MathVector::unit(void)
{
	float s = norm();
	if (s > 0)
	{
		m_vector[0] /= s;
		m_vector[1] /= s;
		m_vector[2] /= s;
	}
	return *this;
}

MathVector MathVector::trunc()
{
	m_vector[0] = floor(m_vector[0]);
	m_vector[1] = floor(m_vector[1]);
	m_vector[2] = floor(m_vector[2]);

	return *this;
}

const float * MathVector::fv(void)
{
	return m_vector;
}

MathVector MathVector::operator +(const MathVector &v)
{
	return MathVector(m_vector[0] + v.m_vector[0], m_vector[1] + v.m_vector[1], m_vector[2] + v.m_vector[2]);
}

MathVector MathVector::operator -(const MathVector &v)
{
	return MathVector(m_vector[0] - v.m_vector[0], m_vector[1] - v.m_vector[1], m_vector[2] - v.m_vector[2]);
}

MathVector MathVector::operator *(const float v)
{
	return MathVector(m_vector[0] * v, m_vector[1] * v, m_vector[2] * v);
}

float MathVector::operator *(const MathVector &v)
{
	return m_vector[0] * v.m_vector[0] + m_vector[1] * v.m_vector[1] + m_vector[2] * v.m_vector[2];
}

MathVector MathVector::operator /(const float v)
{
	return MathVector(m_vector[0] / v, m_vector[1] / v, m_vector[2] / v);
}

MathVector MathVector::operator =(const float *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
	return *this;
}

MathVector MathVector::operator =(const int *v)
{
	m_vector[0] = (float)v[0];
	m_vector[1] = (float)v[1];
	m_vector[2] = (float)v[2];
	return *this;
}

MathVector MathVector::operator =(const float v)
{
	m_vector[0] = v;
	m_vector[1] = v;
	m_vector[2] = v;
	return *this;
}

const float * MathVector::operator ()(void)
{
	return m_vector;
}

float MathVector::operator[] (int id)
{
	return m_vector[id];
}

const MathVector & MathVector::operator +=(const MathVector &v)
{
	m_vector[0] += v.m_vector[0];
	m_vector[1] += v.m_vector[1];
	m_vector[2] += v.m_vector[2];
	return *this;
}

const MathVector & MathVector::operator -=(const MathVector &v)
{
	m_vector[0] -= v.m_vector[0];
	m_vector[1] -= v.m_vector[1];
	m_vector[2] -= v.m_vector[2];
	return *this;
}

const MathVector & MathVector::operator *=(const float v)
{
	m_vector[0] *= v; m_vector[1] *= v; m_vector[2] *= v;
	return *this;
}

const MathVector & MathVector::operator /=(const float v)
{
	m_vector[0] /= v; m_vector[1] /= v; m_vector[2] /= v;
	return *this;
}

bool MathVector::operator ==(const MathVector &v) const
{
	return (m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

bool MathVector::operator !=(const MathVector &v) const
{
	return !(m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

bool MathVector::operator <(const MathVector &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] < v.m_vector[2];
		else
			return m_vector[1] < v.m_vector[1];
	else
		return m_vector[0] < v.m_vector[0];
}

bool MathVector::operator >(const MathVector &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] > v.m_vector[2];
		else
			return m_vector[1] > v.m_vector[1];
	else
		return m_vector[0] > v.m_vector[0];
}

MathVectorD::MathVectorD(void)
{
	m_vector[0] = m_vector[1] = m_vector[2] = 0;
}

MathVectorD::MathVectorD(const double *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
}

MathVectorD::MathVectorD(const double *v1, const double *v2)
{
	m_vector[0] = v2[0] - v1[0];
	m_vector[1] = v2[1] - v1[1];
	m_vector[2] = v2[2] - v1[2];
}

MathVectorD::MathVectorD(const double v1, const double v2, const double v3)
{
	m_vector[0] = v1;
	m_vector[1] = v2;
	m_vector[2] = v3;
}

MathVectorD::MathVectorD(const int v1, const int v2, const int v3)
{
	m_vector[0] = (const double)v1;
	m_vector[1] = (const double)v2;
	m_vector[2] = (const double)v3;
}

MathVectorD::MathVectorD(const MathVectorD &v)
{
	m_vector[0] = v.m_vector[0];
	m_vector[1] = v.m_vector[1];
	m_vector[2] = v.m_vector[2];
}

MathVectorD::~MathVectorD(void)
{
}

MathVectorD MathVectorD::cross(const MathVectorD v)
{
	MathVectorD p;
	p.m_vector[0] = m_vector[1] * v.m_vector[2] - m_vector[2] * v.m_vector[1];
	p.m_vector[1] = m_vector[2] * v.m_vector[0] - m_vector[0] * v.m_vector[2];
	p.m_vector[2] = m_vector[0] * v.m_vector[1] - m_vector[1] * v.m_vector[0];
	return p;
}

double MathVectorD::norm(void)
{
	return sqrt(m_vector[0] * m_vector[0] + m_vector[1] * m_vector[1] + m_vector[2] * m_vector[2]);
}

MathVectorD MathVectorD::unit(void)
{
	double s = norm();
	if (s > 0)
	{
		m_vector[0] /= s;
		m_vector[1] /= s;
		m_vector[2] /= s;
	}
	return *this;
}

MathVectorD MathVectorD::trunc()
{
	m_vector[0] = floor(m_vector[0]);
	m_vector[1] = floor(m_vector[1]);
	m_vector[2] = floor(m_vector[2]);

	return *this;
}

const double * MathVectorD::fv(void)
{
	return m_vector;
}

MathVectorD MathVectorD::operator +(const MathVectorD &v)
{
	return MathVectorD(m_vector[0] + v.m_vector[0], m_vector[1] + v.m_vector[1], m_vector[2] + v.m_vector[2]);
}

MathVectorD MathVectorD::operator -(const MathVectorD &v)
{
	return MathVectorD(m_vector[0] - v.m_vector[0], m_vector[1] - v.m_vector[1], m_vector[2] - v.m_vector[2]);
}

MathVectorD MathVectorD::operator *(const double v)
{
	return MathVectorD(m_vector[0] * v, m_vector[1] * v, m_vector[2] * v);
}

double MathVectorD::operator *(const MathVectorD &v)
{
	return m_vector[0] * v.m_vector[0] + m_vector[1] * v.m_vector[1] + m_vector[2] * v.m_vector[2];
}

MathVectorD MathVectorD::operator /(const double v)
{
	return MathVectorD(m_vector[0] / v, m_vector[1] / v, m_vector[2] / v);
}

MathVectorD MathVectorD::operator =(const double *v)
{
	m_vector[0] = v[0];
	m_vector[1] = v[1];
	m_vector[2] = v[2];
	return *this;
}

MathVectorD MathVectorD::operator =(const int *v)
{
	m_vector[0] = (double)v[0];
	m_vector[1] = (double)v[1];
	m_vector[2] = (double)v[2];
	return *this;
}

MathVectorD MathVectorD::operator =(const double v)
{
	m_vector[0] = v;
	m_vector[1] = v;
	m_vector[2] = v;
	return *this;
}

const double * MathVectorD::operator ()(void)
{
	return m_vector;
}

double MathVectorD::operator[] (int id)
{
	return m_vector[id];
}

const MathVectorD & MathVectorD::operator +=(const MathVectorD &v)
{
	m_vector[0] += v.m_vector[0];
	m_vector[1] += v.m_vector[1];
	m_vector[2] += v.m_vector[2];
	return *this;
}

const MathVectorD & MathVectorD::operator -=(const MathVectorD &v)
{
	m_vector[0] -= v.m_vector[0];
	m_vector[1] -= v.m_vector[1];
	m_vector[2] -= v.m_vector[2];
	return *this;
}

const MathVectorD & MathVectorD::operator *=(const double v)
{
	m_vector[0] *= v; m_vector[1] *= v; m_vector[2] *= v;
	return *this;
}

const MathVectorD & MathVectorD::operator /=(const double v)
{
	m_vector[0] /= v; m_vector[1] /= v; m_vector[2] /= v;
	return *this;
}

bool MathVectorD::operator ==(const MathVectorD &v) const
{
	return (m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

bool MathVectorD::operator !=(const MathVectorD &v) const
{
	return !(m_vector[0] == v.m_vector[0]) && (m_vector[1] == v.m_vector[1]) && (m_vector[2] == v.m_vector[2]);
}

bool MathVectorD::operator <(const MathVectorD &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] < v.m_vector[2];
		else
			return m_vector[1] < v.m_vector[1];
	else
		return m_vector[0] < v.m_vector[0];
}

bool MathVectorD::operator >(const MathVectorD &v) const
{
	if (m_vector[0] == v.m_vector[0])
		if (m_vector[1] == v.m_vector[1])
			return m_vector[2] > v.m_vector[2];
		else
			return m_vector[1] > v.m_vector[1];
	else
		return m_vector[0] > v.m_vector[0];
}

void Coordinate::cart2sph(const double *v, double *phi, double *theta)
{
	// phi: azimuth, theta: elevation
	double d = v[0] * v[0] + v[1] * v[1];
	*phi = (d == 0) ? 0: atan2(v[1], v[0]);
	*theta = (v[2] == 0) ? 0: atan2(v[2], sqrt(d));
}

void Coordinate::cart2sph(const float *v, float *phi, float *theta)
{
	// phi: azimuth, theta: elevation
	float d = v[0] * v[0] + v[1] * v[1];
	*phi = (d == 0) ? 0: atan2(v[1], v[0]);
	*theta = (v[2] == 0) ? 0: atan2(v[2], sqrt(d));
}

void Coordinate::sph2cart(float phi, float theta, float *v)
{
	// phi: azimuth, theta: elevation
	v[2] = sin(theta);
	float coselev = cos(theta);
	v[0] = coselev * cos(phi);
	v[1] = coselev * sin(phi);
}

void Coordinate::sph2cart(double phi, double theta, double *v)
{
	// phi: azimuth, theta: elevation
	v[2] = sin(theta);
	double coselev = cos(theta);
	v[0] = coselev * cos(phi);
	v[1] = coselev * sin(phi);
}

void Coordinate::cart2bary(float *a, float *b, float *c, float *p, float *coeff, float err)
{
	// test dataset for debug
	/*float a[3] = {-0.6498,0.3743,0.6616};
	float b[3] = {-0.6571,0.3837,0.6488};
	float c[3] = {-0.6646,0.3675,0.6506};
	float p[3] = {-0.6572,0.3752,0.6537};
	float coeff[3];*/

	// a counter clockwise order
	MathVector A(a), B(b), C(c), P(p);
	MathVector N((B-A).cross(C-A));

	float ABC = N * N / N.norm();
	N.unit();
	coeff[0] = (B-P).cross(C-P) * N / ABC;
	coeff[1] = (C-P).cross(A-P) * N / ABC;
	//coeff[2] = (A-P).cross(B-P) * N / ABC;
	coeff[2] = 1 - coeff[0] - coeff[1];

	if (fabs(coeff[0]) < err)
	{
		coeff[0] = 0;
		coeff[1] = 1 - (P - B).norm() / (C - B).norm();
		coeff[2] = 1 - coeff[1];
		if (fabs(coeff[1]) < err)
		{
			coeff[1] = 0;
			coeff[2] = 1;
		}
		else if (fabs(coeff[2]) < err)
		{
			coeff[1] = 1;
			coeff[2] = 0;
		}
	}
	else if (fabs(coeff[1]) < err)
	{
		coeff[1] = 0;
		coeff[2] = 1 - (P - C).norm() / (A - C).norm();
		coeff[0] = 1 - coeff[2];
		if (fabs(coeff[2]) < err)
		{
			coeff[2] = 0;
			coeff[0] = 1;
		}
		else if (fabs(coeff[0]) < err)
		{
			coeff[2] = 1;
			coeff[0] = 0;
		}
	}
	else if (fabs(coeff[2]) < err)
	{
		coeff[2] = 0;
		coeff[0] = 1 - (P - A).norm() / (B - A).norm();
		coeff[1] = 1 - coeff[0];
		if (fabs(coeff[0]) < err)
		{
			coeff[0] = 0;
			coeff[1] = 1;
		}
		else if (fabs(coeff[1]) < err)
		{
			coeff[0] = 1;
			coeff[1] = 0;
		}
	}
	// debug
	/*printf("coeff: %f %f %f\n",coeff[0],coeff[1],coeff[2]);
	MathVector PP = A * coeff[0] + B * coeff[1] + C * coeff[2];
	printf("recons: %f %f %f\n", PP.fv()[0],PP.fv()[1],PP.fv()[2]);*/
}
void Coordinate::cart2bary(double *a, double *b, double *c, double *p, double *coeff, double err)
{
	// a counter clockwise order
	MathVectorD A(a), B(b), C(c), P(p);
	MathVectorD N((B-A).cross(C-A));

	double ABC = N * N / N.norm();
	N.unit();
	coeff[0] = (B-P).cross(C-P) * N / ABC;
	coeff[1] = (C-P).cross(A-P) * N / ABC;
	coeff[2] = 1 - coeff[0] - coeff[1];
	
	if (fabs(coeff[0]) < err)
	{
		coeff[0] = 0;
		coeff[1] = 1 - (P - B).norm() / (C - B).norm();
		coeff[2] = 1 - coeff[1];
		if (fabs(coeff[1]) < err)
		{
			coeff[1] = 0;
			coeff[2] = 1;
		}
		else if (fabs(coeff[2]) < err)
		{
			coeff[1] = 1;
			coeff[2] = 0;
		}
	}
	else if (fabs(coeff[1]) < err)
	{
		coeff[1] = 0;
		coeff[2] = 1 - (P - C).norm() / (A - C).norm();
		coeff[0] = 1 - coeff[2];
		if (fabs(coeff[2]) < err)
		{
			coeff[2] = 0;
			coeff[0] = 1;
		}
		else if (fabs(coeff[0]) < err)
		{
			coeff[2] = 1;
			coeff[0] = 0;
		}
	}
	else if (fabs(coeff[2]) < err)
	{
		coeff[2] = 0;
		coeff[0] = 1 - (P - A).norm() / (B - A).norm();
		coeff[1] = 1 - coeff[0];
		if (fabs(coeff[0]) < err)
		{
			coeff[0] = 0;
			coeff[1] = 1;
		}
		else if (fabs(coeff[1]) < err)
		{
			coeff[0] = 1;
			coeff[1] = 0;
		}
	}
}

void Coordinate::rotPoint(const float *p0, const float *mat, float *p1)
{
	for (int i = 0; i < 3; i++)
		p1[i] = mat[i * 3 + 0] * p0[0] + mat[i * 3 + 1] * p0[1] + mat[i * 3 + 2] * p0[2];
}
void Coordinate::rotPointInv(const float *p0, const float *mat, float *p1)
{
	for (int i = 0; i < 3; i++)
		p1[i] = mat[i] * p0[0] + mat[i + 3] * p0[1] + mat[i + 6] * p0[2];
}
void Coordinate::rotPoint(const double *p0, const double *mat, double *p1)
{
	for (int i = 0; i < 3; i++)
		p1[i] = mat[i * 3 + 0] * p0[0] + mat[i * 3 + 1] * p0[1] + mat[i * 3 + 2] * p0[2];
}
void Coordinate::rotPointInv(const double *p0, const double *mat, double *p1)
{
	for (int i = 0; i < 3; i++)
		p1[i] = mat[i] * p0[0] + mat[i + 3] * p0[1] + mat[i + 6] * p0[2];
}
void Coordinate::rotation(const float *axis, const float theta, float *mat)
{
	memset(mat, 0, sizeof(float) * 9);

	mat[0] = 1; mat[4] = 1; mat[8] = 1;

	MathVector RAxis(axis);
	if (RAxis.norm() == 0) return;
	RAxis.unit();

	float A[9] = {0, -RAxis[2], RAxis[1], RAxis[2], 0, -RAxis[0], -RAxis[1], RAxis[0], 0};

	// A * sin(-theta)
	for (int i = 0; i < 9; i++)
		mat[i] += A[i] * sin(theta);
	// A * A * (1 - cos(-theta))
	for (int i = 0; i < 9; i++)
		for (int j = 0; j < 3; j++)
			mat[i] += A[j + (i / 3) * 3] * A[j * 3 + i % 3] * (1 - cos(theta));
}
void Coordinate::rotation(const double *axis, const double theta, double *mat)
{
	memset(mat, 0, sizeof(double) * 9);

	mat[0] = 1; mat[4] = 1; mat[8] = 1;

	MathVectorD RAxis(axis);
	if (RAxis.norm() == 0) return;
	RAxis.unit();

	double A[9] = {0, -RAxis[2], RAxis[1], RAxis[2], 0, -RAxis[0], -RAxis[1], RAxis[0], 0};

	// A * sin(-theta)
	for (int i = 0; i < 9; i++)
		mat[i] += A[i] * sin(theta);
	// A * A * (1 - cos(-theta))
	for (int i = 0; i < 9; i++)
		for (int j = 0; j < 3; j++)
			mat[i] += A[j + (i / 3) * 3] * A[j * 3 + i % 3] * (1 - cos(theta));
}
void Coordinate::rotation2equator(float *v, float *mat, float *pole)
{
	MathVector p0(pole), axis;
	float dot;

	// fit to the equator
	float rot[9];
	MathVector v0(v);
	axis = p0.cross(v0);
	if (axis.norm() == 0) axis = p0;
	dot = p0 * v0;
	dot = (dot > 1) ? 1: dot;
	dot = (dot < -1) ? -1: dot;
	float theta = PI / 2 - acos(dot);
	Coordinate::rotation(axis.fv(), theta, rot);
	float u[3];
	rotPoint(v, rot, u);

	// rotation btw v and u
	MathVector v1(u);
	axis = v0.cross(v1);
	if (axis.norm() == 0) axis = p0;
	dot = v0 * v1;
	dot = (dot > 1) ? 1: dot;
	dot = (dot < -1) ? -1: dot;
	theta = acos(dot);
	rotation(axis.fv(), theta, rot);

	// fit the new pole to the original
	float rot2[9];
	MathVector p1(0, 0, 1);
	axis = p0.cross(p1);
	if (axis.norm() == 0) axis = p0;
	dot = p0 * p1;
	dot = (dot > 1) ? 1: dot;
	dot = (dot < -1) ? -1: dot;
	theta = acos(dot);
	rotation(axis.fv(), theta, rot2);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			MathVector col(rot[j], rot[j + 3], rot[j + 6]);
			MathVector row(rot2[i * 3], rot2[i * 3 + 1], rot2[i * 3 + 2]);
			mat[i * 3 + j] = row * col;
		}
	}

	// debug
	/*printf("Vertex: %f %f %f\n", v[0], v[1], v[2]);
	printf("sph: %f %f\n", phi, theta);
	for (int i = 0; i < 9; i++)
	{
		if (i % 3 == 0) printf("\n");
		printf("%f ", mat[i]);
	}
	printf("\n");*/
}
int Coordinate::intersection(const float *t0, const float *t1, const float *t2, const float a, const float b, const float c, const float d, float *p0, float *p1)
{
	float eps = 1e-8;
	float it0 = a * t0[0] + b * t0[1] + c * t0[2];
	float it1 = a * t1[0] + b * t1[1] + c * t1[2];
	float it2 = a * t2[0] + b * t2[1] + c * t2[2];
	float ait0 = fabs(it0 + d);
	float ait1 = fabs(it1 + d);
	float ait2 = fabs(it2 + d);
	
	// face
	if (ait0 <= eps && ait1 <= eps && ait2 <= eps)
	{
		return 3;
	}
	// vertex
	else if (ait0 <= eps && ait1 > eps && ait2 > eps && (it1 + d) * (it2 + d) > 0)
	{
		p0[0] = p1[0] = t0[0];
		p0[1] = p1[1] = t0[1];
		p0[2] = p1[2] = t0[2];
		return 1;
	}
	else if (ait1 <= eps && ait2 > eps && ait0 > eps && (it2 + d) * (it0 + d) > 0)
	{
		p0[0] = p1[0] = t1[0];
		p0[1] = p1[1] = t1[1];
		p0[2] = p1[2] = t1[2];
		return 1;
	}
	else if (ait2 <= eps && ait0 > eps && ait1 > eps && (it0 + d) * (it1 + d) > 0)
	{
		p0[0] = p1[0] = t2[0];
		p0[1] = p1[1] = t2[1];
		p0[2] = p1[2] = t2[2];
		return 1;
	}
	// edge
	else if (ait0 <= eps && ait1 <= eps && ait2 > eps)
	{
		p0[0] = t0[0]; p0[1] = t0[1]; p0[2] = t0[2];
		p1[0] = t1[0]; p1[1] = t1[1]; p1[2] = t1[2];
		return 2;
	}
	else if (ait1 <= eps && ait2 <= eps && ait0 > eps)
	{
		p0[0] = t1[0]; p0[1] = t1[1]; p0[2] = t1[2];
		p1[0] = t2[0]; p1[1] = t2[1]; p1[2] = t2[2];
		return 2;
	}
	else if (ait2 <= eps && ait0 <= eps && ait1 > eps)
	{
		p0[0] = t2[0]; p0[1] = t2[1]; p0[2] = t2[2];
		p1[0] = t0[0]; p1[1] = t0[1]; p1[2] = t0[2];
		return 2;
	}
	// otherwise
	else
	{
		if ((it0 + d) * (it1 + d) > 0 && (it1 + d) * (it2 + d) > 0 && (it2 + d) * (it0 + d) > 0)
			return -1;

		// intersection at two points
		const float *v0, *v1, *v2;
		float k0, k1, k2;
		float t;
		if ((it1 + d) * (it2 + d) >= 0)
		{
			v0 = t0; v1 = t1; v2 = t2;
			k0 = it0; k1 = it1; k2 = it2;
		}
		else if ((it2 + d) * (it0 + d) >= 0)
		{
			v0 = t1; v1 = t2; v2 = t0;
			k0 = it1; k1 = it2; k2 = it0;
		}
		else if ((it0 + d) * (it1 + d) >= 0)
		{
			v0 = t2; v1 = t0; v2 = t1;
			k0 = it2; k1 = it0; k2 = it1;
		}

		if (fabs(k1 + d) < eps) k1 = -d;
		t = (d + k1) / (k1 - k0);
		p0[0] = v0[0] * t + v1[0] * (1 - t);
		p0[1] = v0[1] * t + v1[1] * (1 - t);
		p0[2] = v0[2] * t + v1[2] * (1 - t);

		if (fabs(k2 + d) < eps) k2 = -d;
		t = (d + k2) / (k2 - k0);
		p1[0] = v0[0] * t + v2[0] * (1 - t);
		p1[1] = v0[1] * t + v2[1] * (1 - t);
		p1[2] = v0[2] * t + v2[2] * (1 - t);

		return 2;
	}
}
void Coordinate::proj2plane(const float a, const float b, const float c, const float d, const float *p0, float *p1)
{
	float portion = (a * p0[0] + b * p0[1] + c * p0[2] + d) / (a * a + b * b + c * c);
	p1[0] = p0[0] - a * portion;
	p1[1] = p0[1] - b * portion;
	p1[2] = p0[2] - c * portion;
}
void Coordinate::proj2plane(const double a, const double b, const double c, const double d, const double *p0, double *p1)
{
	double portion = (a * p0[0] + b * p0[1] + c * p0[2] + d) / (a * a + b * b + c * c);
	p1[0] = p0[0] - a * portion;
	p1[1] = p0[1] - b * portion;
	p1[2] = p0[2] - c * portion;
}
float Coordinate::dpoint2tri(const float *t0, const float *t1, const float *t2, float *p0)
{
	MathVector B(t0), P(p0);
	MathVector E0 = MathVector(t1) - B;
	MathVector E1 = MathVector(t2) - B;

	MathVector D = B - P;
	float a = E0 * E0;
	float b = E0 * E1;
	float c = E1 * E1;
	float d = E0 * D;
	float e = E1 * D;
	float f = D * D;

	float det = a * c - b * b; // do we have to use abs here?
	float s = b * e - c * d;
	float t = b * d - a * e;

	float sqrDistance;

	// Terible tree of conditionals to determine in which region of the diagram
	// shown above the projection of the point into the triangle-plane lies.
	if ((s + t) <= det)
	{
		if (s < 0)
		{
			if (t < 0)
			{
				// region4
				if (d < 0)
				{
					t = 0;
					if (-d >= a)
					{
						s = 1;
						sqrDistance = a + 2 * d + f;
					}
					else
					{
						s = -d / a;
						sqrDistance = d * s + f;
					}
				}
				else
				{
					s = 0;
					if (e >= 0)
					{
						t = 0;
						sqrDistance = f;
					}
					else
					{
						if (-e >= c)
						{
							t = 1;
							sqrDistance = c + 2 * e + f;
						}
						else
						{
							t = -e / c;
							sqrDistance = e * t + f;
						}
					}
				} // end of region 4
			}
			else
			{
				// region 3
				s = 0;
				if (e >= 0)
				{
					t = 0;
					sqrDistance = f;
				}
				else
				{
					if (-e >= c)
					{
						t = 1;
						sqrDistance = c + 2 * e + f;
					}
					else
					{
						t = -e / c;
						sqrDistance = e * t + f;
					}
				}
			} // end of region 3
		}
		else
		{
			if (t < 0)
			{
				// region 5
				t = 0;
				if (d >= 0)
				{
					s = 0;
					sqrDistance = f;
				}
				else
				{
					if (-d >= a)
					{
						s = 1;
						sqrDistance = a + 2 * d + f;
					}
					else
					{
						s = -d / a;
						sqrDistance = d * s + f;
					}
				}
			}
			else
			{
				// region 0
				float invDet = 1 / det;
				s = s * invDet;
				t = t * invDet;
				sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
			}
		}
	}
	else
	{
		if (s < 0)
		{
			// region 2
			float tmp0 = b + d;
			float tmp1 = c + e;
			if (tmp1 > tmp0) // minimum on edge s+t=1
			{
				float numer = tmp1 - tmp0;
				float denom = a - 2 * b + c;
				if (numer >= denom)
				{
					s = 1;
					t = 0;
					sqrDistance = a + 2 * d + f;
				}
				else
				{
					s = numer / denom;
					t = 1 - s;
					sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
				}
			}
			else          // minimum on edge s=0
			{
				s = 0;
				if (tmp1 <= 0)
				{
					t = 1;
					sqrDistance = c + 2 * e + f;
				}
				else
				{
					if (e >= 0)
					{
						t = 0;
						sqrDistance = f;
					}
					else
					{
						t = -e / c;
						sqrDistance = e * t + f;
					}
				}
			} // end of region 2
		}
		else
		{
			if (t < 0)
			{
				// region6 
				float tmp0 = b + e;
				float tmp1 = a + d;
				if (tmp1 > tmp0)
				{
					float numer = tmp1 - tmp0;
					float denom = a - 2 * b + c;
					if (numer >= denom)
					{
						t = 1;
						s = 0;
						sqrDistance = c + 2 * e + f;
					}
					else
					{
						t = numer / denom;
						s = 1 - t;
						sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
					}
				}
				else
				{
					t = 0;
					if (tmp1 <= 0)
					{
						s = 1;
						sqrDistance = a + 2 * d + f;
					}
					else
					{
						if (d >= 0)
						{
							s = 0;
							sqrDistance = f;
						}
						else
						{
							s = -d / a;
							sqrDistance = d * s + f;
						}
					}
				}
			}	// end of region 6
			else
			{
				// region 1
				float numer = c + e - b - d;
				if (numer <= 0)
				{
					s = 0;
					t = 1;
					sqrDistance = c + 2 * e + f;
				}
				else
				{
					float denom = a - 2 * b + c;
					if (numer >= denom)
					{
					  s = 1;
					  t = 0;
					  sqrDistance = a + 2 * d + f;
					}
					else
					{
					  s = numer / denom;
					  t = 1 - s;
					  sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
					}
				} // end of region 1
			}
		}
	}

	// account for numerical round-off error
	if (sqrDistance < 0) sqrDistance = 0;

	float dist = sqrt(sqrDistance);

	return dist;
}
void Coordinate::sphmean(const float *v1, const float *v2, float *v, float w)
{
	if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2])
	{
		v[0] = v1[0]; v[1] = v1[1]; v[2] = v1[2];
	}
	else
	{
		MathVector V1(v1), V2(v2);
		V1.unit(); V2.unit();
		float angle = acos(V1 * V2) * 0.5f * (1.0f - w);
		MathVector axis = V1.cross(V2);
		axis.unit();
		float mat[9];
		rotation(axis.fv(), angle, mat);
		rotPoint(v1, mat, v);

		if (V1.cross(v) * V2.cross(v) > 0)
		{
			rotation(axis.fv(), -angle, mat);
			rotPoint(v1, mat, v);
		}

		MathVector V(v); V.unit();
		v[0] = V[0]; v[1] = V[1]; v[2] = V[2];
	}
}

float Coordinate::arclen(const float *v1, const float *v2)
{
	if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2])
	{
		return 0;
	}
	else
	{
		MathVector V1(v1), V2(v2);
		V1.unit(); V2.unit();

		return fabs(acos(V1 * V2));
	}
}

double Series::factorial(int x)
{
	if (x == 0) return 1;
	return x * factorial(x - 1);
}

double Series::factorial(int x, int stopx)
{
	if (x == stopx) return stopx;
	return x * factorial(x - 1);
}

void Series::legendre(int n, float x, float *Y)
{
	if (n < 0) return;

	float **P = new float*[n + 1];
	for (int i = 0; i <= n; i++) P[i] = new float[n + 1];
	float factor = -sqrt(1.0 - pow(x,2));

	// Init legendre
	P[0][0] = 1.0;        // P_0,0(x) = 1
	if (n == 0)
	{
		Y[0] = P[0][0];
		return;
	}

	// Easy values
	P[1][0] = x;      // P_1,0(x) = x
	P[1][1] = factor;     // P_1,1(x) = −sqrt(1 − x^2)
	if (n == 1)
	{
		Y[0] = P[1][0];
		Y[1] = P[1][1];
		return;
	}

	for (int l = 2; l <= n; l++)
	{
		for (int m = 0; m < l - 1 ; m++)
		{
			// P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
			P[l][m] = ((float)(2 * l - 1) * x * P[l - 1][m] - (float)(l + m - 1) * P[l - 2][m]) / (float)(l - m);
		}
		// P_l,l-1 = (2l-1)*x*P_l-1,l-1
		P[l][l-1] = (float)(2 * l - 1) * x * P[l-1][l-1];
		// P_l,l = (2l-1)*factor*P_l-1,l-1
		P[l][l] = (float)(2 * l - 1) * factor * P[l-1][l-1];
	}

	for (int i = 0; i <= n; i++) Y[i] = P[n][i];

	// release memory
	for (int i = 0; i <= n; i++) delete [] P[i];
	delete [] P;
}

void Series::legendre(int n, double x, double *Y)
{
	if (n < 0) return;

	double **P = new double*[n + 1];
	for (int i = 0; i <= n; i++) P[i] = new double[n + 1];
	double factor = -sqrt(1.0 - pow(x,2));

	// Init legendre
	P[0][0] = 1.0;        // P_0,0(x) = 1
	if (n == 0)
	{
		Y[0] = P[0][0];
		return;
	}

	// Easy values
	P[1][0] = x;      // P_1,0(x) = x
	P[1][1] = factor;     // P_1,1(x) = −sqrt(1 − x^2)
	if (n == 1)
	{
		Y[0] = P[1][0];
		Y[1] = P[1][1];
		return;
	}

	for (int l = 2; l <= n; l++)
	{
		for (int m = 0; m < l - 1 ; m++)
		{
			// P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
			P[l][m] = ((double)(2 * l - 1) * x * P[l - 1][m] - (double)(l + m - 1) * P[l - 2][m]) / (double)(l - m);
		}
		// P_l,l-1 = (2l-1)*x*P_l-1,l-1
		P[l][l-1] = (double)(2 * l - 1) * x * P[l-1][l-1];
		// P_l,l = (2l-1)*factor*P_l-1,l-1
		P[l][l] = (double)(2 * l - 1) * factor * P[l-1][l-1];
	}

	for (int i = 0; i <= n; i++) Y[i] = P[n][i];

	// release memory
	for (int i = 0; i <= n; i++) delete [] P[i];
	delete [] P;
}

float Statistics::sum(float *A, int n)
{
	float res = 0;
	for (int i = 0; i < n; i++) res += A[i];
	return res;
}
double Statistics::sum(double *A, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++) res += A[i];
	return res;
}
float Statistics::mean(float *A, int n)
{
	return sum(A, n) / (float)n;
}
double Statistics::mean(double *A, int n)
{
	return sum(A, n) / (double)n;
}

float Statistics::min(float *A, int n)
{
	float res = A[0];
	for (int i = 1; i < n; i++)
		if (res > A[i]) res = A[i];
	return res;
}
double Statistics::min(double *A, int n)
{
	double res = A[0];
	for (int i = 1; i < n; i++)
		if (res > A[i]) res = A[i];
	return res;
}

float Statistics::max(float *A, int n)
{
	float res = A[0];
	for (int i = 1; i < n; i++)
		if (res < A[i]) res = A[i];
	return res;
}
double Statistics::max(double *A, int n)
{
	double res = A[0];
	for (int i = 1; i < n; i++)
		if (res < A[i]) res = A[i];
	return res;
}
void Statistics::sum(float *A, int n, int dim, float *res)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < dim; j++)
			res[j] += A[i * dim + j];
}
void Statistics::sum(double *A, int n, int dim, double *res)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < dim; j++)
			res[j] += A[i * dim + j];
}
void Statistics::mean(float *A, int n, int dim, float *res)
{
	sum(A, n, dim, res);
	for (int i = 0; i < dim; i++)
		res[i] /= (float)n;
}
void Statistics::mean(double *A, int n, int dim, double *res)
{
	sum(A, n, dim, res);
	for (int i = 0; i < dim; i++)
		res[i] /= (double)n;
}
float Statistics::correlation(float *A, float *B, int n)
{
	float res = 0;
	float meanA = mean(A, n);
	float meanB = mean(B, n);

	for (int i = 0; i < n; i++)
		res += (A[i] - meanA) * (B[i] - meanB);

	return res;
}
float Statistics::var(float *A, int n)
{
	float m = mean(A, n);
	float res = 0;
	for (int i = 0; i < n; i++) res += (A[i] - m) * (A[i] - m) / n;
	return res;
}
double Statistics::var(double *A, int n)
{
	double m = mean(A, n);
	double res = 0;
	for (int i = 0; i < n; i++) res += (A[i] - m) * (A[i] - m) / n;
	return res;
}
void Statistics::cov(float *p, int n, int dim, float *M)
{
	float *x;
	float *m = new float[dim];

	memset(m, 0, sizeof(float) * dim);

	// mean
	mean(p, n, dim, m);

	// covariance matrix (spd)
	for (int i = 0; i < n; i++)
	{
		x = &p[i * dim];
		for (int j = 0; j < dim; j++)
		{
			for (int k = j; k < dim; k++)
			{
				M[j * dim + k] += (x[j] - m[j]) * (x[k] - m[k]) / (float)(n - 1);
				M[k * dim + j] = M[j * dim + k];
			}
		}
	}

	delete m;
}
void Statistics::cov_trans(float *p, int n, int dim, float *M)
{
	float *x, *y;
	float *m = new float[dim];

	memset(m, 0, sizeof(float) * dim);

	// mean
	mean(p, n, dim, m);

	// covariance matrix (spd)
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			x = &p[i * dim];
			y = &p[j * dim];
			for (int k = 0; k < dim; k++)
			{
				M[i * n + j] += (x[k] - m[k]) * (y[k] - m[k]) / (float)(n - 1);
			}
			M[j * n + i] = M[i * n + j];
		}
	}

	delete m;
}
void Statistics::wcov_trans(float *p, int n, int dim, float *M, float *w)
{
	float *x, *y;
	float *m = new float[dim];

	memset(m, 0, sizeof(float) * dim);

	// mean
	mean(p, n, dim, m);

	// covariance matrix (spd)
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			x = &p[i * dim];
			y = &p[j * dim];
			for (int k = 0; k < dim; k++)
			{
				M[i * n + j] += (x[k] - m[k]) * (y[k] - m[k]) * w[k] * w[k] / (float)(n - 1);
			}
			M[j * n + i] = M[i * n + j];
		}
	}

	delete m;
}
void Statistics::wcov_trans(double *p, int n, int dim, double *M, double *w)
{
	double *x, *y;
	double *m = new double[dim];

	memset(m, 0, sizeof(double) * dim);

	// mean
	mean(p, n, dim, m);

	// covariance matrix (spd)
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			x = &p[i * dim];
			y = &p[j * dim];
			for (int k = 0; k < dim; k++)
			{
				M[i * n + j] += (x[k] - m[k]) * (y[k] - m[k]) * w[k] * w[k] / (double)(n - 1);
			}
			M[j * n + i] = M[i * n + j];
		}
	}

	delete m;
}
float Statistics::NCC(float *A, float *B, int n)
{
	return correlation(A, B, n) / sqrt(correlation(A, A, n)) / sqrt(correlation(B, B, n));
}
float Statistics::normal_pdf(float x, float u, float sigma)
{
	return exp( -1 * (x - u) * (x - u) / (2 * sigma * sigma)) / (sigma * sqrt(2 * PI));
}
float Statistics::normal_cdf_approx(float x, float u, float sigma)
{
	const double ninf = u - 10 * sigma;
	double sum = 0;
	double n = 1e3; // tune for speed/accuracy
	double c = (x - ninf) / n;

	for (int k = 1; k < n - 1; k++)
		sum += normal_pdf(ninf + (float)k * c, u, sigma);

	return c * ((normal_pdf(x, u, sigma) + normal_pdf(ninf, u, sigma)) / 2 + sum);
}
float Statistics::erf(float x)
{
	float sign = (x > 0) ? 1: -1;
	x *= sign;
	float y = 1.0 / ( 1.0 + 0.3275911 * x);
	float e = 1 - (((((
			+ 1.061405429  * y
			- 1.453152027) * y
			+ 1.421413741) * y
			- 0.284496736) * y
			+ 0.254829592) * y)
			* exp (-x * x);
	return e * sign;
}
float Statistics::normal_cdf(float x, float u, float sigma)
{
	return 0.5 * (1 + erf((x - u) / (sigma * sqrt(2.))));
}
float Statistics::median(float *v, int n)
{
	// copy to tmp
	float *tmp = new float[n];
	memcpy(tmp, v, sizeof(float) * n);
	// sort
	std::sort(tmp, tmp + n);
	//median
	float med;
	if (n % 2 == 1) med = tmp[n / 2];
	else med = (tmp[n / 2 - 1] + tmp[n / 2]) * 0.5f;

	delete [] tmp;

	return med;
}
double Statistics::median(double *v, int n)
{
	// copy to tmp
	double *tmp = new double[n];
	memcpy(tmp, v, sizeof(double) * n);
	// sort
	std::sort(tmp, tmp + n);
	//median
	float med;
	if (n % 2 == 1) med = tmp[n / 2];
	else med = (tmp[n / 2 - 1] + tmp[n / 2]) * 0.5f;

	delete [] tmp;

	return med;
}

void LinearAlgebra::eig3symmetric(const double M[3][3], double lambda[3], double eigv[3][3])
{
	// http://en.wikipedia.org/wiki/Eigenvalue_algorithm
	double p1 = M[0][1] * M[0][1] + M[0][2] * M[0][2] + M[1][2] * M[1][2];
	if (p1 == 0)	// diagonal matrix
	{
		if (M[0][0] < M[1][1] && M[0][0] < M[2][2])
		{
			lambda[0] = M[0][0];
			if (M[1][1] < M[2][2])
			{
				lambda[1] = M[1][1];
				lambda[2] = M[2][2];
				eigv[0][0] = 1; eigv[0][1] = 0; eigv[0][2] = 0;
				eigv[1][0] = 0; eigv[1][1] = 1; eigv[1][2] = 0;
				eigv[2][0] = 0; eigv[2][1] = 0; eigv[2][2] = 1;
			}
			else
			{
				lambda[1] = M[2][2];
				lambda[2] = M[1][1];
				eigv[0][0] = 1; eigv[0][1] = 0; eigv[0][2] = 0;
				eigv[1][0] = 0; eigv[1][1] = 0; eigv[1][2] = 1;
				eigv[2][0] = 0; eigv[2][1] = 1; eigv[2][2] = 0;
			}
		}
		else if (M[1][1] < M[2][2] && M[1][1] < M[0][0])
		{
			lambda[0] = M[1][1];
			if (M[0][0] < M[2][2])
			{
				lambda[1] = M[0][0];
				lambda[2] = M[2][2];
				eigv[0][0] = 0; eigv[0][1] = 1; eigv[0][2] = 0;
				eigv[1][0] = 1; eigv[1][1] = 0; eigv[1][2] = 0;
				eigv[2][0] = 0; eigv[2][1] = 0; eigv[2][2] = 1;
			}
			else
			{
				lambda[1] = M[2][2];
				lambda[2] = M[0][0];
				eigv[0][0] = 0; eigv[0][1] = 1; eigv[0][2] = 0;
				eigv[1][0] = 0; eigv[1][1] = 0; eigv[1][2] = 1;
				eigv[2][0] = 1; eigv[2][1] = 0; eigv[2][2] = 0;
			}
		}
		else if (M[2][2] < M[0][0] && M[2][2] < M[1][1])
		{
			lambda[0] = M[2][2];
			if (M[0][0] < M[1][1])
			{
				lambda[1] = M[0][0];
				lambda[2] = M[1][1];
				eigv[0][0] = 0; eigv[0][1] = 0; eigv[0][2] = 1;
				eigv[1][0] = 1; eigv[1][1] = 0; eigv[1][2] = 0;
				eigv[2][0] = 0; eigv[2][1] = 1; eigv[2][2] = 0;
			}
			else
			{
				lambda[1] = M[1][1];
				lambda[2] = M[0][0];
				eigv[0][0] = 0; eigv[0][1] = 0; eigv[0][2] = 1;
				eigv[1][0] = 0; eigv[1][1] = 1; eigv[1][2] = 0;
				eigv[2][0] = 1; eigv[2][1] = 0; eigv[2][2] = 0;
			}
		}
	}
	else
	{
		double q = trace3(M) / 3; // trace(M)
		double p2 = (M[0][0] - q) * (M[0][0] - q) + (M[1][1] - q) * (M[1][1] - q) + (M[2][2] - q) * (M[2][2] - q) + 2 * p1;
		double p = sqrt(p2 / 6);
		double A[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				A[i][j] = (1 / p) * (M[i][j] - q * (int)(i == j));
		double r = det3(A) / 2;

		// In exact arithmetic for a symmetric matrix  -1 <= r <= 1
		// but computation error can leave it slightly outside this range.
		double phi;
		if (r <= -1) phi = PI / 3;
		else if (r >= 1) phi = 0;
		else phi = acos(r) / 3;

		// the eigenvalues satisfy eig3 <= eig2 <= eig1
		lambda[0] = q + 2 * p * cos(phi + (2 * PI / 3));	// the smallest eigenvalue
		lambda[2] = q + 2 * p * cos(phi);
		lambda[1] = 3 * q - lambda[0] - lambda[2];	// since trace(A) = eig1 + eig2 + eig3

		double **B = new double*[3];
		double *work = new double[9];
		for (int i = 0; i < 3; i++) B[i] = &work[i * 3];
		// M - lambda * I
		for (int n = 0; n < 3; n++)
		{
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					B[i][j] = M[i][j] - lambda[n] * (int)(i == j);
			echelonEigv((const double **)B, 3, eigv[n]);
			double norm = 0;
			for (int i = 0; i < 3; i++) norm += eigv[n][i] * eigv[n][i];
			norm = sqrt(norm);
			for (int i = 0; i < 3; i++) eigv[n][i] /= norm; 
		}
		
		delete [] B;
		delete [] work;
	}
}
double LinearAlgebra::trace3(const double M[3][3])
{
	return M[0][0] + M[1][1] + M[2][2];
}
double LinearAlgebra::det3(const double M[3][3])
{
	return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - 
			M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + 
			M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}

void LinearAlgebra::echelonEigv(const double **M, int n, double *x)
{
	// square matrix only
	double **A = new double*[n];	// rows
	double *work = new double[n * n];
	for (int i = 0; i < n; i++)
	{
		A[i] = &work[n * i];
		for (int j = 0; j < n; j++)
			A[i][j] = M[i][j];
	}
	
	for (int i = 0; i < n - 1; i++)	// row
	{
		for (int j = i + 1; j < n; j++) // lower rows
		{
			double c = M[j][i] / M[i][i];
			for (int k = 0; k < n; k++)
				A[j][k] -= A[i][k] * c;
		}
	}

	// solution
	x[n - 1] = 1;
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = 0;
		for (int j = i + 1; j < n; j++)
			x[i] += A[i][j] * x[j];
		x[i] /= -A[i][i];
	}

	delete [] A;
	delete [] work;
}

