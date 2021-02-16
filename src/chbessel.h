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
#include "qgaus.h"
#include <cmath>
#include <cstdint>
#include <iostream>



namespace FastBessel {

#define CALC_CHEB_SUM(x, sum) \
	{						  \
		double z = d/x; \
		double z2 = 2.*(2*z-1.); \
		double ch_0 = 1.; \
		double ch_1 = (2.*z-1.); \
		double ch_2; \
		sum = ch_0*ck[0] + ch_1*ck[1]; \
		for (int i = 2; i <= M; ++i) { \
			ch_2 = z2*ch_1 - ch_0; \
			sum += ch_2*ck[i]; \
			ch_0 = ch_1; \
			ch_1 = ch_2; \
		} \
	};

#define POWER_SUM(x, sum) \
	{ \
		long double xO2 = 0.5*x; \
		const long double A = -(log(xO2) + EUL_GAMMA_LD); \
		long double xO2sq = xO2*xO2; \
		for (int i = MAXIT_IKBESS-1; i>= 1; i--) { \
			sum = xO2sq*(coef[i]*(A+1./(2.*i+1.)+ns[i]) + sum); \
		} \
		sum = x*A+x+x*sum; \
	}

#define FAST_CHEB_SUM(x, sum) \
		{ \
		double z = d/x; \
		sum = 1.25331413731550e+00+ \
			z*(-3.91660667890412e-01+ \
			z*(3.15776408289850e-01+ \
			z*(-4.06194424195295e-01+ \
			z*(7.19597732003498e-01+ \
			z*(-1.62708622337386e+00+ \
			z*(4.46581915059184e+00+ \
			z*(-1.42188765206382e+01+ \
			z*(4.97110608920295e+01+ \
			z*(-1.78790459190384e+02+ \
			z*(6.21256867237806e+02+ \
			z*(-1.98596426994130e+03+ \
			z*(5.65054653532254e+03+ \
			z*(-1.40148864771337e+04+ \
			z*(2.99068912288121e+04+ \
			z*(-5.44194703279499e+04+ \
			z*(8.38397009322609e+04+ \
			z*(-1.08611516820013e+05+ \
			z*(1.17395578367076e+05+ \
			z*(-1.04843244388150e+05+ \
			z*(7.63708826316075e+04+ \
			z*(-4.45748586206830e+04+ \
			z*(2.03233348896703e+04+ \
			z*(-6.96738114384029e+03+ \
			z*(1.68804707483831e+03+ \
			z*(-2.57529748166367e+02+ \
			z*1.85968591391546e+01)))))))))))))))))))))))));\
		};

static const double PI = 3.141592653589793;
static const long double PI2 = 1.5707'9632'6794'8966'2;
static const long double PI2D = 1.5707'9632'6794'897;
static const double EUL_GAMMA_D = 0.5772'1566'4901'5329;
static const long double EUL_GAMMA_LD = 0.5772'1566'4901'5328'6;
static const __float128 PI2Q = 1.5707963267948966192313216916397514q;
static const __float128 SQRT_PI2Q = 1.2533141373155002512078826424055226q;
static const double TINY = std::numeric_limits<double>::min();

class Bess {
	// calculates modified bessel function of second kind and its integral using Chebyshev polynomial expansion
public:
	Bess(const bool fast = false, const int n = 34, const int m = 34, const double dd = 2.);
	//-----------
	double ik00x_ch(const double x) const; // t: [0, x], chebyshev approximation for x >= d;
	long double ik00x_pwr(const double x) const; // t: [0, x], power approximation for t < d;
	double ik00x(const double x) const; // adaptive for any x > 0
	//-----------
	double ik0ab_ch(const double x1, const double x2) const; // t: [a,b], chebyshev approximation for a,b >= d;
	long double ik0ab_pwr(const double x1, const double x2) const; // t: [a,b], power approximation for a,b < d;
	double ik0ab_num(const double x1, const double x2) const; // numerical solution using gauss quadrature
	double ik0ab(const double x1, const double x2) const; // adaptive function uses _ch, _pwr and _num functions;
	double abs_ik0ab(const double x1, const double x2) const;
	//-----------
	double _k0(const double x) const; // sqrt(x)*exp(x)*BesselK[0,x], x >= d;
	double k0(const double x) const; // BesselK[0,x], x >= d
private:
	const int MAXIT_IKBESS=20;
	const int GAUSS_POINTS = 20;
	const bool is_fast;
	const double d;
	const int N, M;
	std::vector<long double>  _coef(), _ns();
	std::vector<double> fak(const int m, const int n, double d); // calculates coefficients of Chebyshev series for _k0(x)
	std::vector<double> fck(const int m, const int n, double d, double mu); // calculates coefficients of Chebyshev series for _ik0(x)
	std::vector<__float128> fakq(const int m, const int n, const __float128 d, const __float128 multt = 1.q); // same as fak but with quadruple presicion
	std::vector<__float128> fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt = 1.q); // same as fck but with quadruple precision
	const std::vector<double> ak, ck; // coefficients for constructor Bess()
	const std::vector<long double> coef, ns;
	//BessK0 bess;
	//GaussIntegrator gs;
	std::vector<double> xgs, wgs; // gauss abscissas and weights
	static const double k0pi[5],k0qi[3],k0p[5],k0q[3],k0pp[8],k0qq[8];
	inline double poly(const double* cof, const int n, const double x) const; // Evaluate a polynomial

};

}
