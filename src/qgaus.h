/*
 * qgaus.h
 *
 *  Created on: 26 но€б. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <cmath>

void gauleg(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w);

struct GaussIntegrator {
	int n;
	std::vector<double> x, w;
	explicit GaussIntegrator(const int n_);
	template <class T>
	double Integrate(const T& func, const double a, const double b) const {
		double xm=0.5*(b+a);
		double xr=0.5*(b-a);
		double s=0;
		int jmax = x.size()/2;
		for (int j=0;j<jmax;j++) {
			double dx=xr*x[j];
			s += w[j]*(func(xm+dx)+func(xm-dx));
		}
		return s *= xr;
	};
};

template <class T>
double qgaus(T& func, const double a, const double b)
//Returns the integral of the function or functor func between a and b, by ten-point Gauss-
//Legendre integration: the function is evaluated exactly ten times at interior points in the range
//of integration.
{
//Here are the abscissas and weights:
	static const double x[]={0.1488743389816312,0.4333953941292472,
							 0.6794095682990244,0.8650633666889845,0.9739065285171717};
	static const double w[]={0.2955242247147529,0.2692667193099963,
								0.2190863625159821,0.1494513491505806,0.0666713443086881};
	double xm=0.5*(b+a);
	double xr=0.5*(b-a);
	double s=0; //Will be twice the average value of the function, since the
				//ten weights (five numbers above each used twice)
				//sum to 2.
	for (int j=0;j<5;j++) {
		double dx=xr*x[j];
		s += w[j]*(func(xm+dx)+func(xm-dx));
	}
	return s *= xr; //Scale the answer to the range of integration.
}
