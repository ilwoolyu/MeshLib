/*************************************************
*	SphericalHarmonics.cpp
*
*	Release: August 2012
*	Update: June 2021
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*	
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
*************************************************/

#include <cmath>
#include "SphericalHarmonics.h"

void SphericalHarmonics::basis(int degree, const float *p, float *Y, int from)
{
	// real spherical harmonics basis functions
	// polar coordinate
	double phi, theta;
	double dp[3] = {p[0], p[1], p[2]};
	Coordinate::cart2sph(dp, &phi, &theta);
	theta = PI / 2 - theta;  // convert to interval [0, PI]
	double *Pm = new double[degree + 1];

	// square root of 2
	const double sqr2 = sqrt(2.0f);

	double *P[2];
	P[0] = new double[degree + 1];
	P[1] = new double[degree + 1];
	int base = 0;

	for (int l = from; l <= degree; l++)
	{
		// legendre part
		bool schmidt = l > 85;
		if (!schmidt)
		{
			Series::legendre(l, cos(theta), Pm, P, base);
			base = l + 1;
		}
		else
		{
			Series::legendre2(l, cos(theta), Pm, schmidt);
		}
		double lconstant = sqrt((2 * l + 1) / (4 * PI));

		int center = (l + 1) * (l + 1) - l - 1 - from * from;

		Y[center] = (float)(lconstant * Pm[0]);

		for (int m = 1; m <= l; m++)
		{
			double precoeff = (schmidt) ? lconstant: sqr2 * lconstant * sqrt(1.0 / Series::factorial(l + m, l - m + 1));

			if (!schmidt && m % 2 == 1) precoeff = -precoeff;
			Y[center + m] = (float)(precoeff * Pm[m] * cos(m * phi));
			Y[center - m] = (float)(precoeff * Pm[m] * sin(m * phi));
		}
	}

	delete [] P[0];
	delete [] P[1];
	delete [] Pm;
}

void SphericalHarmonics::basis(int degree, const double *p, double *Y, int from)
{
	// real spherical harmonics basis functions
	// polar coordinate
	double phi, theta;
	Coordinate::cart2sph(p, &phi, &theta);
	theta = PI / 2 - theta;  // convert to interval [0, PI]
	double *Pm = new double[degree + 1];

	// square root of 2
	const double sqr2 = sqrt(2.0);

	double *P[2];
	P[0] = new double[degree + 1];
	P[1] = new double[degree + 1];
	int base = 0;

	for (int l = from; l <= degree; l++)
	{
		// legendre part
		bool schmidt = l > 85;
		if (!schmidt)
		{
			Series::legendre(l, cos(theta), Pm, P, base);
			base = l + 1;
		}
		else
		{
			Series::legendre2(l, cos(theta), Pm, schmidt);
		}
		double lconstant = sqrt((2 * l + 1) / (4 * PI));

		int center = (l + 1) * (l + 1) - l - 1 - from * from;

		Y[center] = lconstant * Pm[0];

		for (int m = 1; m <= l; m++)
		{
			double precoeff = (schmidt) ? lconstant: sqr2 * lconstant * sqrt(1.0 / Series::factorial(l + m, l - m + 1));

			if (!schmidt && m % 2 == 1) precoeff = -precoeff;
			Y[center + m] = precoeff * Pm[m] * cos(m * phi);
			Y[center - m] = precoeff * Pm[m] * sin(m * phi);
		}
	}

	delete [] P[0];
	delete [] P[1];
	delete [] Pm;
}
