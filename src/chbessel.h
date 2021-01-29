/*
 * chbessel.h
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <quadmath.h>
#include <limits>
#include "bessel.h"
#include "qgaus.h"
#include <cmath>

namespace FastBessel {

static const double PI = 3.141592653589793;
static const long double PI2 = 1.5707'9632'6794'8966'2;
static const double EUL_GAMMA_D = 0.5772'1566'4901'5329;
static const long double EUL_GAMMA_LD = 0.5772'1566'4901'5328'6;
static const __float128 PI2Q = 1.5707963267948966192313216916397514q;
static const __float128 SQRT_PI2Q = 1.2533141373155002512078826424055226q;
static const double TINY = std::numeric_limits<double>::min();


struct Bess {
	// calculates modified bessel function of second kind and its integral using Chebyshev polynomial expansion
	Bess(const int n = 34, const int m = 34, const double dd = 2.);
	const int MAXIT_IKBESS=20;
	const int GAUSS_POINTS = 20;
	const double d;
	const int N, M;
	double _k0(const double x); // sqrt(x)*exp(x)*BesselK[0,x], x >= d;
	double _ik0(const double x) const; // sqrt(x)*exp(x)*Integral[BesselK[x, inf]], x >= d;
	double k0(const double x); // BesselK[0,x], x >= d
	double ik0(const double x) const; // Integral BessK0[x, inf], x >= d
	double ik02(const double x) const; // Integral BessK0[x, inf], x < 2
	long double _ik02(const double x) const; // Integral BessK0[0, x], x < 2

	double analytic_ik0ab(const double x1, const double x2) const; // integral BessK0[x1, x2]
	double ik0ab(const double x1, const double x2) const;
	std::vector<long double>  _coef(), _ns();
	std::vector<double> fak(const int m, const int n, double d); // calculates coefficients of Chebyshev series for _k0(x)
	std::vector<double> fck(const int m, const int n, double d, double mu); // calculates coefficients of Chebyshev series for _ik0(x)
	std::vector<__float128> fakq(const int m, const int n, const __float128 d, const __float128 multt = 1.q); // same as fak but with quadruple presicion
	std::vector<__float128> fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt = 1.q); // same as fck but with quadruple precision
	const std::vector<double> ak, ck; // coefficients for constructor Bess()
	std::vector<double> cha, chc;
	const std::vector<long double> coef, ns;
	BessK0 bess;
	GaussIntegrator gs;
	double num_ik0ab(const double x1, const double x2) const;
};

}
