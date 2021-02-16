/*
 * quadrature.h
 *
 *  Created on: 26 íîÿá. 2020 ã.
 *      Author: Dmitry_Di
 */

#pragma once

#include "interp_1d.h"
#include <cmath>
#include <vector>

struct Quadrature {
// Abstract base class for elementary quadrature algorithms
	int n;
	virtual double next() = 0;
};

template <class T>
struct Trapzd: Quadrature {
	// Routine implementing the extended trapezoidal rule.
	T& func;
	double a, b, s;

	Trapzd(T& funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb), s(0) {n=0;}

	double next() {
		double x,tnm,sum,del;
		int it,j;
		n++;
		if (n == 1) {
			s=0.5*(b-a)*(func(a)+func(b));
			return s;
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
};

template <class T>
double qtrap(T& func, const double a, const double b, const double eps=1.0e-10) {
	//Returns the integral of the function or functor func from a to b. The constants EPS can be
	//set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
	//allowed number of steps. Integration is performed by the trapezoidal rule.
	const int JMAX=100;
	double s,olds=0.0; 						//Initial value of olds is arbitrary.
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		s=t.next();
		if (j > 5) 	{						//Avoid spurious early convergence.
			if (std::abs(s-olds) < eps*std::abs(olds) ||
					(s == 0.0 && olds == 0.0))
				return s;
		}
		olds=s;
	}
	throw("Too many steps in routine qtrap");
}

template<class T>
double qsimp(T &func, const double a, const double b, const double eps=1.0e-10) {
//Returns the integral of the function or functor func from a to b. The constants EPS can be
//set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
//allowed number of steps. Integration is performed by Simpson’s rule.
	const int JMAX=20;
	double s,st,ost=0.0,os=0.0;
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		st=t.next();
		s=(4.0*st-ost)/3.0;
		if (j > 5) { //Avoid spurious early convergence.
			if (std::abs(s-os) < eps*std::abs(os) ||
					(s == 0.0 && os == 0.0))
				return s;
		}
		os=s;
		ost=st;
	}
	throw("Too many steps in routine qsimp");
}

template <class T>
double qromb(T &func, double a, double b, const double eps=1.0e-10) {
//Returns the integral of the function or functor func from a to b. Integration is performed by
//Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.
	const int JMAX=20, JMAXP=JMAX+1, K=5;
	//Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
	//JMAX limits the total number of steps; K is the number of points used in the
	//extrapolation.
	std::vector<double> s(JMAX),h(JMAXP); //These store the successive trapezoidal approxi
	Poly_interp polint(h,s,K);     //mations and their relative stepsizes.
	h[0]=1.0;
	Trapzd<T> t(func,a,b);
	for (int j=1;j<=JMAX;j++) {
		s[j-1]=t.next();
		if (j >= K) {
			double ss=polint.rawinterp(j-K,0.0);
			if (std::abs(polint.dy) <= eps*std::abs(ss)) return ss;
		}
		h[j]=0.25*h[j-1];
		//This is a key step: The factor is 0.25 even though the stepsize is decreased by only
		//0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1),
		//not just a polynomial in h.
	}
	throw std::runtime_error("Too many steps in routine qromb");
}


template <class T>
struct Midpnt : Quadrature {
	//Routine implementing the extended midpoint rule.
	double a,b,s; //Limits of integration and current value of integral
	T &funk;
	Midpnt(T& funcc, const double aa, const double bb) :
	funk(funcc), a(aa), b(bb) {n=0;}
//The constructor takes as inputs func, the function or functor to be integrated between
//limits a and b, also input.
	double next(){
	//Returns the nth stage of refinement of the extended midpoint rule. On the first call (n=1),
	//the routine returns the crudest estimate of
	//R b
	//a f.x/dx. Subsequent calls set n=2,3,... and
	//improve the accuracy by adding .2=3/  3n-1 additional interior points.
		int it,j;
		double x,tnm,sum,del,ddel;
		n++;
		if (n == 1) {
			return (s=(b-a)*func(0.5*(a+b)));
		} else {
			for(it=1,j=1;j<n-1;j++) it *= 3;
				tnm=it;
				del=(b-a)/(3.0*tnm);
				ddel=del+del; //The added points alternate in spacing be
				x=a+0.5*del;  //tween del and ddel.
				sum=0.0;
				for (j=0;j<it;j++) {
					sum += func(x);
					x += ddel;
					sum += func(x);
					x += del;
				}
				s=(s+(b-a)*sum/tnm)/3.0; //The new sum is combined with the old inte
				return s; 				 //gral to give a refined integral.
		}
	}
	virtual double func(const double x) {return funk(x);} //Identity mapping.
};

template<class T>
double qromo(Midpnt<T>& q, const double eps=3.0e-9) {
//Romberg integration on an open interval. Returns the integral of a function using any specified
//elementary quadrature algorithm q and Romberg’s method. Normally q will be an open formula,
//not evaluating the function at the endpoints. It is assumed that q triples the number of steps
//on each call, and that its error series contains only even powers of the number of steps. The
//routines midpnt, midinf, midsql, midsqu, midexp are possible choices for q. The constants
//below have the same meanings as in qromb.
	const int JMAX=14, JMAXP=JMAX+1, K=5;
	std::vector<double> h(JMAXP),s(JMAX);
	Poly_interp polint(h,s,K);
	h[0]=1.0;
	for (int j=1;j<=JMAX;j++) {
		s[j-1]=q.next();
		if (j >= K) {
			double ss=polint.rawinterp(j-K,0.0);
			if (std::abs(polint.dy) <= eps*std::abs(ss)) return ss;
		}
		h[j]=h[j-1]/9.0; //This is where the assumption of step tripling and an even
	} 				//error series is used.
	throw std::runtime_error("Too many steps in routine qromo");
}
