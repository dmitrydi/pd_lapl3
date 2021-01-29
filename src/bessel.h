/*
 * bessel.h
 *
 *  Created on: 25 но€б. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <quadmath.h>
#include <utility>

namespace FastBessel {

struct Bessik {
	Bessik();
	static const double i0p[14],i0q[5],i0pp[5],i0qq[6];
	static const double i1p[14],i1q[5],i1pp[5],i1qq[6];
	static const double k0pi[5],k0qi[3],k0p[5],k0q[3],k0pp[8],k0qq[8];
	static const double k1pi[5],k1qi[3],k1p[5],k1q[3],k1pp[8],k1qq[8];
	double y,z,ax,term;

	double i0(const double x);	// Returns the modified Bessel function I0(x) for any real x

	double i1(const double x);	// Returns the modified Bessel function I1(x) for any real x

	double k0(const double x) const;	// Returns the modified Bessel function K0(x) for positive real x.


	double k1(const double x); // Returns the modified Bessel function K1(x) for positive real x.

	double in(const int n, const double x);
	// Returns the modified Bessel function In(x) for any real x and n >= 0;

	double kn(const int n, const double x);
	//Returns the modified Bessel function Kn(x) for positive x and n >= 0

	inline double poly(const double* cof, const int n, const double x) const; // Evaluate a polynomial
};

class BessK0 {
private:
	Bessik bess;
public:
	BessK0();
	double operator()(const double x) const;
};
}
