/*************************************************
*	Geom.h
*
*	Release: July 2011
*	Update: June 2021
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
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
	MathVector operator -(void);
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
	MathVectorD operator -(void);
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
	static void cart2bary(const float *a, const float *b, const float *c, const float *p, float *coeff, float err = 0);
	static void cart2bary(const double *a, const double *b, const double *c, const double *p, double *coeff, double err = 0);
	static void rotPoint(const float *p0, const float *mat, float *p1);
	static void rotPoint(const double *p0, const double *mat, double *p1);
	static void rotPointInv(const float *p0, const float *mat, float *p1);
	static void rotPointInv(const double *p0, const double *mat, double *p1);
	static void rotation(const float *axis, const float theta, float *mat);
	static void rotation(const double *axis, const double theta, double *mat);
	static void rotation2equator(const float *v, float *mat, const float *pole);
	static int intersection(const float *t0, const float *t1, const float *t2, const float a, const float b, const float c, const float d, float *p0, float *p1);
	static void proj2plane(const float a, const float b, const float c, const float d, const float *p0, float *p1);
	static void proj2plane(const double a, const double b, const double c, const double d, const double *p0, double *p1);
	static float dpoint2tri(const float *t0, const float *t1, const float *t2, const float *p0);
	static float dpoint2tri(const float *t0, const float *t1, const float *t2, const float *p0, float *coeff);
	static void sphmean(const float *v1, const float *v2, float *v, float w = 0.5f);
	static float arclen(const float *v1, const float *v2);
};

class Series
{
public:
	static double factorial(int x);
	static double factorial(int x, int stopx);
	static void legendre(int n, float x, float *Y, float **preP = NULL, int base = 0);
	static void legendre(int n, double x, double *Y, double **preP = NULL, int base = 0);
	static void legendre2(int n, float x, float *Y, bool schmidt = false);
	static void legendre2(int n, double x, double *Y, bool schmidt = false);
};

class Statistics
{
public:
	static float sum(const float *A, int n);
	static double sum(const double *A, int n);
	static float mean(const float *A, int n);
	static double mean(const double *A, int n);
	static float min(const float *A, int n);
	static double min(const double *A, int n);
	static float max(const float *A, int n);
	static double max(const double *A, int n);
	static void sum(const float *A, int n, int dim, float *res);
	static void sum(const double *A, int n, int dim, double *res);
	static void mean(const float *A, int n, int dim, float *res);
	static void mean(const double *A, int n, int dim, double *res);
	static float correlation(const float *A, const float *B, int n);
	static float var(const float *A, int n);
	static double var(const double *A, int n);
	static void cov(const float *p, int n, int dim, float *M);
	static void cov_trans(const float *p, int n, int dim, float *M);
	static void wcov_trans(const float *p, int n, int dim, float *M, const float *w);
	static void wcov_trans(const double *p, int n, int dim, double *M, const double *w);
	static float NCC(const float *A, const float *B, int n);
	static float normal_pdf(float x, float u, float sigma);
	static float normal_cdf_approx(float x, float u, float sigma);
	static float erf(float x);
	static float normal_cdf(float x, float u, float sigma);
	static float median(const float *v, int n);
	static double median(const double *v, int n);
};

class LinearAlgebra
{
public:
	static void eig3symmetric(const double M[3][3], double lambda[3], double eigv[3][3]);
	static void Ab(const double **A, int n, int m, const double *b, double *x);
	static double trace3(const double M[3][3]);
	static double det3(const double M[3][3]);
private:
	static void gaussElim(const double **M, int n, double *x);
	static void echelonEigv(const double **M, int n, double *x);
};

#endif
