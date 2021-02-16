/*
 * interp_1d.h
 *
 *  Created on: 26 но€б. 2020 г.
 *      Author: Dmitry_Di
 *   	Contains: utilities for function interpolation
 */
#pragma once

#include <cmath>
#include <vector>

struct Base_interp {
	//Abstract class used for interpolation
	int n, mm, jsav, cor, dj;
	const double *xx, *yy;
	Base_interp(const std::vector<double>& x,
			const double *y, int m);
	double interp(double x);	// Given a value x, return an interpolated value, using data pointed to by xx and yy.
	int locate(const double x);
	int hunt(const double x);
	double virtual rawinterp(int jlo, double x) = 0;
	//Derived classes provide this as the actual interpolation method.
};

struct Linear_interp : Base_interp
/*
 Piecewise linear interpolation object. Construct with x and y vectors, then call interp for
 interpolated values.

 USE:
 int n=...;
 vector<double> xx(n), yy(n);
 ...
 Linear_interp myfunc(xx,yy);
 double x,y;
 ...
 y = myfunc.interp(x);

*/
{
	Linear_interp(const std::vector<double>& xv, const std::vector<double>& yv);
	double rawinterp(int j, double x);
};

struct Poly_interp : Base_interp
/*
Polynomial interpolation object. Construct with x and y vectors, and the number M of points
to be used locally (polynomial order plus one), then call interp for interpolated values.
*/
{
	double dy;
	Poly_interp(const std::vector<double> &xv, const std::vector<double> &yv, int m);
	double rawinterp(int jl, double x);
};
