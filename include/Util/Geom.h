/*************************************************
*	Geom.h
*
*	Release: July 2011
*	Update: October 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#ifndef GEOM_HH_
#define GEOM_HH_

#define PI 3.141592653589793238462643

class MathVector
{
public:
	MathVector(void);
	MathVector(const float *v);
	MathVector(const float *v1, const float *v2);
	MathVector(const float v1, const float v2, const float v3);
	MathVector(const int v1, const int v2, const int v3);
	MathVector(const MathVector &v);
	~MathVector(void);
	MathVector cross(const MathVector v);
	float norm(void);
	MathVector unit(void);
	MathVector trunc();
	const float * fv(void);
	MathVector operator +(const MathVector &v);
	MathVector operator -(const MathVector &v);
	MathVector operator *(const float v);
	float operator *(const MathVector &v);
	MathVector operator /(const float v);
	MathVector operator =(const float *v);
	MathVector operator =(const int *v);
	MathVector operator =(const float v);
	const float *operator ()(void);
	float operator[] (int id);
	const MathVector & operator +=(const MathVector &v);
	const MathVector & operator -=(const MathVector &v);
	const MathVector & operator *=(const float v);
	const MathVector & operator /=(const float v);
	bool operator ==(const MathVector &v) const;
	bool operator !=(const MathVector &v) const;
	bool operator <(const MathVector &v) const;
	bool operator >(const MathVector &v) const;
private:
	float m_vector[3];
};

class MathVectorD
{
public:
	MathVectorD(void);
	MathVectorD(const double *v);
	MathVectorD(const double *v1, const double *v2);
	MathVectorD(const double v1, const double v2, const double v3);
	MathVectorD(const int v1, const int v2, const int v3);
	MathVectorD(const MathVectorD &v);
	~MathVectorD(void);
	MathVectorD cross(const MathVectorD v);
	double norm(void);
	MathVectorD unit(void);
	MathVectorD trunc();
	const double * fv(void);
	MathVectorD operator +(const MathVectorD &v);
	MathVectorD operator -(const MathVectorD &v);
	MathVectorD operator *(const double v);
	double operator *(const MathVectorD &v);
	MathVectorD operator /(const double v);
	MathVectorD operator =(const double *v);
	MathVectorD operator =(const int *v);
	MathVectorD operator =(const double v);
	const double *operator ()(void);
	double operator[] (int id);
	const MathVectorD & operator +=(const MathVectorD &v);
	const MathVectorD & operator -=(const MathVectorD &v);
	const MathVectorD & operator *=(const double v);
	const MathVectorD & operator /=(const double v);
	bool operator ==(const MathVectorD &v) const;
	bool operator !=(const MathVectorD &v) const;
	bool operator <(const MathVectorD &v) const;
	bool operator >(const MathVectorD &v) const;
private:
	double m_vector[3];
};


class Coordinate
{
public:
	static void cart2sph(const double *v, double *phi, double *theta);
	static void cart2sph(const float *v, float *phi, float *theta);
	static void sph2cart(float phi, float theta, float *v);
	static void sph2cart(double phi, double theta, double *v);
	static void cart2bary(float *a, float *b, float *c, float *p, float *coeff, float err = 0);
	static void cart2bary(double *a, double *b, double *c, double *p, double *coeff, double err = 0);
	static void rotPoint(const float *p0, const float *mat, float *p1);
	static void rotPoint(const double *p0, const double *mat, double *p1);
	static void rotPointInv(const float *p0, const float *mat, float *p1);
	static void rotPointInv(const double *p0, const double *mat, double *p1);
	static void rotation(const float *axis, const float theta, float *mat);
	static void rotation(const double *axis, const double theta, double *mat);
	static void rotation2equator(float *v, float *mat, float *pole);
	static int intersection(const float *t0, const float *t1, const float *t2, const float a, const float b, const float c, const float d, float *p0, float *p1);
	static void proj2plane(const float a, const float b, const float c, const float d, const float *p0, float *p1);
	static void proj2plane(const double a, const double b, const double c, const double d, const double *p0, double *p1);
	static float dpoint2tri(const float *t0, const float *t1, const float *t2, float *p0);
	static void sphmean(const float *v1, const float *v2, float *v, float w = 0.5f);
	static float arclen(const float *v1, const float *v2);
};

class Series
{
public:
	static double factorial(int x);
	static double factorial(int x, int stopx);
	static void legendre(int n, float x, float *Y);
	static void legendre(int n, double x, double *Y);
};

class Statistics
{
public:
	static float sum(float *A, int n);
	static double sum(double *A, int n);
	static float mean(float *A, int n);
	static double mean(double *A, int n);
	static float min(float *A, int n);
	static double min(double *A, int n);
	static float max(float *A, int n);
	static double max(double *A, int n);
	static void sum(float *A, int n, int dim, float *res);
	static void sum(double *A, int n, int dim, double *res);
	static void mean(float *A, int n, int dim, float *res);
	static void mean(double *A, int n, int dim, double *res);
	static float correlation(float *A, float *B, int n);
	static float var(float *A, int n);
	static double var(double *A, int n);
	static void cov(float *p, int n, int dim, float *M);
	static void cov_trans(float *p, int n, int dim, float *M);
	static void wcov_trans(float *p, int n, int dim, float *M, float *w);
	static void wcov_trans(double *p, int n, int dim, double *M, double *w);
	static float NCC(float *A, float *B, int n);
	static float normal_pdf(float x, float u, float sigma);
	static float normal_cdf_approx(float x, float u, float sigma);
	static float erf(float x);
	static float normal_cdf(float x, float u, float sigma);
	static float median(float *v, int n);
	static double median(double *v, int n);
};

class LinearAlgebra
{
public:
	static void eig3symmetric(const double M[3][3], double lambda[3], double eigv[3][3]);
	static double trace3(const double M[3][3]);
	static double det3(const double M[3][3]);
private:
	static void echelonEigv(const double **M, int n, double *x);
};

#endif
