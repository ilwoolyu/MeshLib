/*************************************************
*	SphericalHarmonics.cpp
*
*	Release: August 2012
*	Update: September 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <cmath>
#include "SphericalHarmonics.h"

void SphericalHarmonics::basis(int degree, float *p, float *Y, int from)
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

	for (int l = from; l <= degree; l++)
	{
		// legendre part
		Series::legendre(l, cos(theta), Pm);
		double lconstant = sqrt((2 * l + 1) / (4 * PI));

		int center = (l + 1) * (l + 1) - l - 1 - from * from;

		Y[center] = (float)(lconstant * Pm[0]);

		for (int m = 1; m <= l; m++)
		{
			double precoeff = lconstant * sqrt(1.0 / Series::factorial(l + m, l - m + 1));

			if (m % 2 == 1) precoeff = -precoeff;
			Y[center + m] = (float)(sqr2 * precoeff * Pm[m] * cos(m * phi));
			Y[center - m] = (float)(sqr2 * precoeff * Pm[m] * sin(m * phi));
		}
	}

	delete [] Pm;
}

void SphericalHarmonics::basis(int degree, double *p, double *Y, int from)
{
	// real spherical harmonics basis functions
	// polar coordinate
	double phi, theta;
	Coordinate::cart2sph(p, &phi, &theta);
	theta = PI / 2 - theta;  // convert to interval [0, PI]
	double *Pm = new double[degree + 1];

	// square root of 2
	const double sqr2 = sqrt(2.0);

	for (int l = from; l <= degree; l++)
	{
		// legendre part
		Series::legendre(l, cos(theta), Pm);
		double lconstant = sqrt((2 * l + 1) / (4 * PI));

		int center = (l + 1) * (l + 1) - l - 1 - from * from;

		Y[center] = lconstant * Pm[0];

		for (int m = 1; m <= l; m++)
		{
			double precoeff = lconstant * sqrt(1.0 / Series::factorial(l + m, l - m + 1));

			if (m % 2 == 1) precoeff = -precoeff;
			Y[center + m] = sqr2 * precoeff * Pm[m] * cos(m * phi);
			Y[center - m] = sqr2 * precoeff * Pm[m] * sin(m * phi);
		}
	}

	delete [] Pm;
}

