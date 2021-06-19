/*************************************************
*	Geom.cpp
*
*	Release: July 2011
*	Update: June 2021
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*	
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
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

MathVector MathVector::operator -(void)
{
	return MathVector(-m_vector[0], -m_vector[1], -m_vector[2]);
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

MathVectorD MathVectorD::operator -(void)
{
	return MathVectorD(-m_vector[0], -m_vector[1], -m_vector[2]);
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

void Coordinate::cart2bary(const float *a, const float *b, const float *c, const float *p, float *coeff, float err)
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
void Coordinate::cart2bary(const double *a, const double *b, const double *c, const double *p, double *coeff, double err)
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
void Coordinate::rotation2equator(const float *v, float *mat, const float *pole)
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
float Coordinate::dpoint2tri(const float *t0, const float *t1, const float *t2, const float *p0, float *coeff)
{
	// projection to the triangle
	float p0_proj[3];
	MathVector N = MathVector(t0, t1).cross(MathVector(t1, t2)).unit();
	Coordinate::proj2plane(N[0], N[1], N[2], -N * MathVector(t0), p0, p0_proj);

	// bary centric
	Coordinate::cart2bary(t0, t1, t2, p0_proj, coeff);

	// contacting point
	MathVector P;
	float eps = 0;
	if (coeff[0] < eps || coeff[1] < eps || coeff[2] < eps)
	{
		float len_e_ab = MathVector(t0, t1).norm();
		float len_e_bc = MathVector(t1, t2).norm();
		float len_e_ca = MathVector(t2, t0).norm();

		float len_ab = MathVector(t0, t1).unit() * MathVector(t0, p0);
		float len_bc = MathVector(t1, t2).unit() * MathVector(t1, p0);
		float len_ca = MathVector(t2, t0).unit() * MathVector(t2, p0);

		if (len_ab > len_e_ab) len_ab = len_e_ab;
		if (len_bc > len_e_bc) len_bc = len_e_bc;
		if (len_ca > len_e_ca) len_ca = len_e_ca;
		if (len_ab < 0) len_ab = 0;
		if (len_bc < 0) len_bc = 0;
		if (len_ca < 0) len_ca = 0;

		MathVector V_ab = MathVector(t0, t1).unit() * len_ab + MathVector(t0);
		MathVector V_bc = MathVector(t1, t2).unit() * len_bc + MathVector(t1);
		MathVector V_ca = MathVector(t2, t0).unit() * len_ca + MathVector(t2);

		float dist_ab = (V_ab - p0).norm();
		float dist_bc = (V_bc - p0).norm();
		float dist_ca = (V_ca - p0).norm();

		if (dist_ab <= dist_bc && dist_ab <= dist_ca)
		{
			coeff[1] = len_ab / len_e_ab;
			coeff[0] = 1 - coeff[1];
			coeff[2] = 0;
			P = V_ab;
		}
		else if (dist_bc <= dist_ca && dist_bc <= dist_ab)
		{
			coeff[2] = len_bc / len_e_bc;
			coeff[1] = 1 - coeff[2];
			coeff[0] = 0;
			P = V_bc;
		}
		else
		{
			coeff[0] = len_ca / len_e_ca;
			coeff[2] = 1 - coeff[0];
			coeff[1] = 0;
			P = V_ca;
		}
	}
	else
	{
		P = MathVector(p0_proj);
	}
	return (MathVector(p0) - P).norm();
}
float Coordinate::dpoint2tri(const float *t0, const float *t1, const float *t2, const float *p0)
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
	const double pre[16] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
							39916800, 479001600, 6227020800, 87178291200, 1307674368000};
	if (x <= 15) return pre[x];
	double fac = pre[15];
	for (int i = 16; i <= x; i++) fac *= i;
	return fac;
}

double Series::factorial(int x, int stopx)
{
	const double pre[16] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
							39916800, 479001600, 6227020800, 87178291200, 1307674368000};
	if (x <= 15 && stopx <= 15) return pre[x] / pre[stopx - 1];
	double fac = 1;
	for (int i = stopx; i <= x; i++) fac *= i;
	return fac;
}

void Series::legendre(int n, float x, float *Y, float **preP, int base)
{
	if (n < 0) return;

	if (n == 0)
	{
		Y[0] = 1.0;
		if (preP != NULL)
			preP[0][0] = 1.0;
		return;
	}

	float factor = -sqrt(1.0 - x * x);
	if (n == 1)
	{
		Y[0] = x;
		Y[1] = factor;
		if (preP != NULL)
		{
			preP[0][0] = 1.0;
			preP[1][0] = x;
			preP[1][1] = factor;
		}
		return;
	}

	float **P;
	bool pre = (preP != NULL && base > 1);
	if (!pre)
	{
		base = 2;
		P = new float*[n + 1];
		for (int i = 0; i <= n; i++) P[i] = new float[i + 1];
		P[0][0] = 1.0;		// P_0,0(x) = 1
		P[1][0] = x;		// P_1,0(x) = x
		P[1][1] = factor;	// P_1,1(x) = −sqrt(1 − x^2)
	}
	else
	{
		P = new float*[n + 1 - base + 2];
		for (int i = 0; i <= n - base + 2; i++) P[i] = new float[i + 1 + base - 2];
		memcpy(P[0], preP[0], sizeof(float) * (base - 1));
		memcpy(P[1], preP[1], sizeof(float) * base);
	}

	for (int l = base; l <= n; l++)
	{
		for (int m = 0; m < l - 1 ; m++)
		{
			// P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
			P[l - base + 2][m] = ((float)(2 * l - 1) * x * P[l - 1 - base + 2][m] - (float)(l + m - 1) * P[l - 2 - base + 2][m]) / (float)(l - m);
		}
		// P_l,l-1 = (2l-1)*x*P_l-1,l-1
		P[l - base + 2][l - 1] = (float)(2 * l - 1) * x * P[l - 1 - base + 2][l - 1];
		// P_l,l = (2l-1)*factor*P_l-1,l-1
		P[l - base + 2][l] = (float)(2 * l - 1) * factor * P[l - 1 - base + 2][l - 1];
	}

	for (int i = 0; i <= n; i++) Y[i] = P[n - base + 2][i];

	// release memory
	if (preP != NULL)
	{
		if (pre)
			std::swap(preP[0], preP[1]);
		else
			memcpy(preP[0], P[n - 1 - base + 2], sizeof(float) * n);
		memcpy(preP[1], P[n - base + 2], sizeof(float) * (n + 1));
	}
	for (int i = 0; i <= n - base + 2; i++) delete [] P[i];
	delete [] P;
}

void Series::legendre(int n, double x, double *Y, double **preP, int base)
{
	if (n < 0) return;

	if (n == 0)
	{
		Y[0] = 1.0;
		if (preP != NULL)
			preP[0][0] = 1.0;
		return;
	}

	double factor = -sqrt(1.0 - x * x);
	if (n == 1)
	{
		Y[0] = x;
		Y[1] = factor;
		if (preP != NULL)
		{
			preP[0][0] = 1.0;
			preP[1][0] = x;
			preP[1][1] = factor;
		}
		return;
	}

	double **P;
	bool pre = (preP != NULL && base > 1);
	if (!pre)
	{
		base = 2;
		P = new double*[n + 1];
		for (int i = 0; i <= n; i++) P[i] = new double[i + 1];
		P[0][0] = 1.0;		// P_0,0(x) = 1
		P[1][0] = x;		// P_1,0(x) = x
		P[1][1] = factor;	// P_1,1(x) = −sqrt(1 − x^2)
	}
	else
	{
		P = new double*[n + 1 - base + 2];
		for (int i = 0; i <= n - base + 2; i++) P[i] = new double[i + 1 + base - 2];
		memcpy(P[0], preP[0], sizeof(double) * (base - 1));
		memcpy(P[1], preP[1], sizeof(double) * base);
	}

	for (int l = base; l <= n; l++)
	{
		for (int m = 0; m < l - 1 ; m++)
		{
			// P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
			P[l - base + 2][m] = ((double)(2 * l - 1) * x * P[l - 1 - base + 2][m] - (double)(l + m - 1) * P[l - 2 - base + 2][m]) / (double)(l - m);
		}
		// P_l,l-1 = (2l-1)*x*P_l-1,l-1
		P[l - base + 2][l - 1] = (double)(2 * l - 1) * x * P[l - 1 - base + 2][l - 1];
		// P_l,l = (2l-1)*factor*P_l-1,l-1
		P[l - base + 2][l] = (double)(2 * l - 1) * factor * P[l - 1 - base + 2][l - 1];
	}

	for (int i = 0; i <= n; i++) Y[i] = P[n - base + 2][i];

	// release memory
	if (preP != NULL)
	{
		if (pre)
			std::swap(preP[0], preP[1]);
		else
			memcpy(preP[0], P[n - 1 - base + 2], sizeof(double) * n);
		memcpy(preP[1], P[n - base + 2], sizeof(double) * (n + 1));
	}
	for (int i = 0; i <= n - base + 2; i++) delete [] P[i];
	delete [] P;
}

void Series::legendre2(int n, float x, float *Y, bool schmidt)
{
	// source from https://www.mathworks.com/help/matlab/ref/legendre.html
	if (n < 0) return;

	if (n == 0)
	{
		Y[0] = 1;
		return;
	}

	float factor = sqrt(1.0 - x * x);
	if (n == 1)
	{
		Y[0] = x;
		Y[1] = (schmidt) ? factor: -factor;
		return;
	}

	float *rootn = new float[2 * n + 1];
	for (int i = 0; i < 2 * n + 1; i++) rootn[i] = sqrt((float)i);
	float twocot = (factor == 0) ? 0: -2 * x / factor;
	float sn = pow(-factor, n);
	memset(Y, 0, sizeof(float) * (n + 1));

	const float tol = 1.0842022e-19;
	if (factor > 0 && fabs(sn) <= tol)
	{
		// Approx solution of x*ln(x) = y
		float v = 9.2 - log(tol) / (n * factor);
		float w = 1 / log(v);
		float m1 = 1 + n * factor * v * w * (1.0058 + w * (3.819 - w * 12.173));

		int mm1 = std::min(n, (int)floor(m1));

		// Start recursion with proper sign
		const float tstart = 1.1920929e-07;
		Y[mm1 - 1] = (x < 0) ? (n + 1) % 2: mm1 % 2;
		Y[mm1 - 1] = (Y[mm1 - 1] == 0) ? -tstart: tstart;

		// Recur from m1 to m = 0, accumulating normalizing factor.
		float sumsq = tol;
		for (int m = mm1 - 2; m >= 0; m--)
		{
			Y[m] = ((m + 1) * twocot * Y[m + 1] - rootn[n + m + 2] * rootn[n - m - 1] * Y[m + 2]) / (rootn[n + m + 1] * rootn[n - m]);
			sumsq += Y[m] * Y[m];
		}
		float scale = 1 / sqrt(2.0f * sumsq - Y[0] * Y[0]);
		for (int i = 0; i <= mm1; i++)
			Y[i] *= scale;
	}

	if (x != 1 && fabs(sn) >= tol)
	{
		float c = 1;
		for (int i = 1; i < n + 1; i++)
			c *= (1.0 - 1.0 / (i * 2));
		Y[n] = sqrt(c) * sn;
		Y[n - 1] = Y[n] * twocot * n / rootn[2 * n];

		for (int m = n - 2; m >= 0; m--)
			Y[m] = (Y[m + 1] * twocot * (m + 1) - Y[m + 2] * rootn[n + m + 2] * rootn[n - m - 1]) / (rootn[n + m + 1] * rootn[n - m]);
	}

	if (factor == 0)
		Y[0] = pow(x, n);

	if (!schmidt)
	{
		for (int m = 1; m < n; m++)
			for (int j = n - m + 1; j < n + m + 1; j++)
				Y[m] *= rootn[j];
		for (int j = 1; j < 2 * n + 1; j++)
			Y[n] *= rootn[j];
	}
	else
	{
		float const1 = -1;
		for (int j = 1; j < n + 1; j++)
		{
			Y[j] *= rootn[2] * const1;
			const1 *= -1;
		}
	}
	delete [] rootn;
}

void Series::legendre2(int n, double x, double *Y, bool schmidt)
{
	// source from https://www.mathworks.com/help/matlab/ref/legendre.html
	if (n < 0) return;

	if (n == 0)
	{
		Y[0] = 1;
		return;
	}

	double factor = sqrt(1.0 - x * x);
	if (n == 1)
	{
		Y[0] = x;
		Y[1] = (schmidt) ? factor: -factor;
		return;
	}

	double *rootn = new double[2 * n + 1];
	for (int i = 0; i < 2 * n + 1; i++) rootn[i] = sqrt((double)i);
	double twocot = (factor == 0) ? 0: -2 * x / factor;
	double sn = pow(-factor, n);
	memset(Y, 0, sizeof(double) * (n + 1));

	const double tol = 1.491668146240041e-154;
	if (factor > 0 && fabs(sn) <= tol)
	{
		// Approx solution of x*ln(x) = y
		double v = 9.2 - log(tol) / (n * factor);
		double w = 1 / log(v);
		double m1 = 1 + n * factor * v * w * (1.0058 + w * (3.819 - w * 12.173));

		int mm1 = std::min(n, (int)floor(m1));

		// Start recursion with proper sign
		const double tstart = 2.22044604925031308e-16;
		Y[mm1 - 1] = (x < 0) ? (n + 1) % 2: mm1 % 2;
		Y[mm1 - 1] = (Y[mm1 - 1] == 0) ? -tstart: tstart;

		// Recur from m1 to m = 0, accumulating normalizing factor.
		double sumsq = tol;
		for (int m = mm1 - 2; m >= 0; m--)
		{
			Y[m] = ((m + 1) * twocot * Y[m + 1] - rootn[n + m + 2] * rootn[n - m - 1] * Y[m + 2]) / (rootn[n + m + 1] * rootn[n - m]);
			sumsq += Y[m] * Y[m];
		}
		double scale = 1 / sqrt(2.0 * sumsq - Y[0] * Y[0]);
		for (int i = 0; i <= mm1; i++)
			Y[i] *= scale;
	}

	if (x != 1 && fabs(sn) >= tol)
	{
		double c = 1;
		for (int i = 1; i < n + 1; i++)
			c *= (1.0 - 1.0 / (i * 2));
		Y[n] = sqrt(c) * sn;
		Y[n - 1] = Y[n] * twocot * n / rootn[2 * n];

		for (int m = n - 2; m >= 0; m--)
			Y[m] = (Y[m + 1] * twocot * (m + 1) - Y[m + 2] * rootn[n + m + 2] * rootn[n - m - 1]) / (rootn[n + m + 1] * rootn[n - m]);
	}

	if (factor == 0)
		Y[0] = pow(x, n);

	if (!schmidt)
	{
		for (int m = 1; m < n; m++)
			for (int j = n - m + 1; j < n + m + 1; j++)
				Y[m] *= rootn[j];
		for (int j = 1; j < 2 * n + 1; j++)
			Y[n] *= rootn[j];
	}
	else
	{
		double const1 = -1;
		for (int j = 1; j < n + 1; j++)
		{
			Y[j] *= rootn[2] * const1;
			const1 *= -1;
		}
	}
	delete [] rootn;
}

float Statistics::sum(const float *A, int n)
{
	float res = 0;
	for (int i = 0; i < n; i++) res += A[i];
	return res;
}
double Statistics::sum(const double *A, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++) res += A[i];
	return res;
}
float Statistics::mean(const float *A, int n)
{
	return sum(A, n) / (float)n;
}
double Statistics::mean(const double *A, int n)
{
	return sum(A, n) / (double)n;
}

float Statistics::min(const float *A, int n)
{
	float res = A[0];
	for (int i = 1; i < n; i++)
		if (res > A[i]) res = A[i];
	return res;
}
double Statistics::min(const double *A, int n)
{
	double res = A[0];
	for (int i = 1; i < n; i++)
		if (res > A[i]) res = A[i];
	return res;
}

float Statistics::max(const float *A, int n)
{
	float res = A[0];
	for (int i = 1; i < n; i++)
		if (res < A[i]) res = A[i];
	return res;
}
double Statistics::max(const double *A, int n)
{
	double res = A[0];
	for (int i = 1; i < n; i++)
		if (res < A[i]) res = A[i];
	return res;
}
void Statistics::sum(const float *A, int n, int dim, float *res)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < dim; j++)
			res[j] += A[i * dim + j];
}
void Statistics::sum(const double *A, int n, int dim, double *res)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < dim; j++)
			res[j] += A[i * dim + j];
}
void Statistics::mean(const float *A, int n, int dim, float *res)
{
	sum(A, n, dim, res);
	for (int i = 0; i < dim; i++)
		res[i] /= (float)n;
}
void Statistics::mean(const double *A, int n, int dim, double *res)
{
	sum(A, n, dim, res);
	for (int i = 0; i < dim; i++)
		res[i] /= (double)n;
}
float Statistics::correlation(const float *A, const float *B, int n)
{
	float res = 0;
	float meanA = mean(A, n);
	float meanB = mean(B, n);

	for (int i = 0; i < n; i++)
		res += (A[i] - meanA) * (B[i] - meanB);

	return res;
}
float Statistics::var(const float *A, int n)
{
	float m = mean(A, n);
	float res = 0;
	for (int i = 0; i < n; i++) res += (A[i] - m) * (A[i] - m) / n;
	return res;
}
double Statistics::var(const double *A, int n)
{
	double m = mean(A, n);
	double res = 0;
	for (int i = 0; i < n; i++) res += (A[i] - m) * (A[i] - m) / n;
	return res;
}
void Statistics::cov(const float *p, int n, int dim, float *M)
{
	const float *x;
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
void Statistics::cov_trans(const float *p, int n, int dim, float *M)
{
	const float *x, *y;
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
void Statistics::wcov_trans(const float *p, int n, int dim, float *M, const float *w)
{
	const float *x, *y;
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
void Statistics::wcov_trans(const double *p, int n, int dim, double *M, const double *w)
{
	const double *x, *y;
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
float Statistics::NCC(const float *A, const float *B, int n)
{
	return correlation(A, B, n) / sqrt(correlation(A, A, n)) / sqrt(correlation(B, B, n));
}
float Statistics::normal_pdf(const float x, const float u, float sigma)
{
	return exp( -1 * (x - u) * (x - u) / (2 * sigma * sigma)) / (sigma * sqrt(2 * PI));
}
float Statistics::normal_cdf_approx(const float x, const float u, float sigma)
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
float Statistics::median(const float *v, int n)
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
double Statistics::median(const double *v, int n)
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
void LinearAlgebra::Ab(const double **A, int n, int m, const double *b, double *x)
{
	// least squares: ATA needs to be full rank to gaurantee a unique solution
	// A = n by m, n >= m
	// b = n by 1
	double *work = new double[m * m];
	double **ATA = new double*[m];
	double *ATb = x;

	for (int i = 0; i < m; i++) ATA[i] = &work[m * i];
	for (int i = 0; i < m; i++)
	{
		for (int j = i; j < m; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				ATA[i][j] += A[k][i] * A[k][j];
			}
			ATA[j][i] = ATA[i][j];
		}
		ATb[i] = 0;
		for (int j = 0; j < n; j++)
		{
			ATb[i] += A[j][i] * b[j];
		}
	}
	gaussElim((const double **)ATA, m, ATb);

	delete [] work;
	delete [] ATA;
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
void LinearAlgebra::gaussElim(const double **M, int n, double *x)
{
	// square matrix only
	double **A = new double*[n];	// rows
	double *work = new double[n * (n + 1)];
	double *tRow = new double[n + 1];
	for (int i = 0; i < n; i++)
	{
		A[i] = &work[(n + 1) * i];
		for (int j = 0; j < n; j++)
			A[i][j] = M[i][j];
		A[i][n] = x[i];
	}

	for (int i = 0; i < n - 1; i++)	// row
	{
		// maximum pivot
		double maxElem = fabs(A[i][i]);
		int r = i;
		for (int j = i + 1; j < n; j++) // lower rows
		{
			double absA = fabs(A[j][i]);
			if (maxElem < absA)
			{
				maxElem = absA;
				r = j;
			}
		}
		// swap
		if (r != i)
		{
			memcpy(tRow, A[r], sizeof(double) * (n + 1));
			memcpy(A[r], A[i], sizeof(double) * (n + 1));
			memcpy(A[i], tRow, sizeof(double) * (n + 1));
			double tx = x[r];
			x[r] = x[i];
			x[i] = tx;
		}
		// elimination
		for (int j = i + 1; j < n; j++) // lower rows
		{
			double c = A[j][i] / A[i][i];
			A[j][i] = 0;
			for (int k = i + 1; k <= n; k++)
				A[j][k] -= A[i][k] * c;
		}
	}

	// solution
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = A[i][n];
		for (int j = i + 1; j < n; j++)
			x[i] -= A[i][j] * x[j];
		x[i] /= A[i][i];
	}

	delete [] A;
	delete [] work;
	delete [] tRow;
}
void LinearAlgebra::echelonEigv(const double **M, int n, double *x)
{
	// square matrix only
	double **A = new double*[n];	// rows
	double *work = new double[n * n];
	double *tRow = new double[n];
	for (int i = 0; i < n; i++)
	{
		A[i] = &work[n * i];
		for (int j = 0; j < n; j++)
			A[i][j] = M[i][j];
	}
	
	for (int i = 0; i < n - 2; i++)	// row
	{
		// maximum pivot
		double maxElem = fabs(A[i][i]);
		int r = i;
		for (int j = i + 1; j < n; j++) // lower rows
		{
			double absA = fabs(A[j][i]);
			if (maxElem < absA)
			{
				maxElem = absA;
				r = j;
			}
		}
		// swap
		if (r != i)
		{
			memcpy(tRow, A[r], sizeof(double) * n);
			memcpy(A[r], A[i], sizeof(double) * n);
			memcpy(A[i], tRow, sizeof(double) * n);
		}
		// elimination
		for (int j = i + 1; j < n - 1; j++) // lower rows
		{
			double c = A[j][i] / A[i][i];
			A[j][i] = 0;
			for (int k = i + 1; k < n; k++)
				A[j][k] -= A[i][k] * c;
		}
	}

	// solution
	x[n - 1] = 1;	// fix the last elem to 1
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = 0;
		for (int j = i + 1; j < n; j++)
			x[i] -= A[i][j] * x[j];
		x[i] /= A[i][i];
	}

	delete [] A;
	delete [] work;
	delete [] tRow;
}

