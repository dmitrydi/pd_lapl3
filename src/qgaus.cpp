/*
 * qgaus.cpp
 *
 *  Created on: 21 ÿíâ. 2021 ã.
 *      Author: Dmitry_Di
 */

#include "qgaus.h"

using namespace std;


GaussIntegrator::GaussIntegrator(const int n_): n(n_), x(n), w(n) {
	gauleg(-1., 1., x, w);
};

void gauleg(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w)
//Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
//and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point
//quadrature formula.
{
	const double EPS=1.0e-14; //EPS is the relative precision.
	double z1,z,xm,xl,pp,p3,p2,p1;
	int n=x.size();
	int m=(n+1)/2; //The roots are symmetric in the interval, so
	xm=0.5*(x2+x1); //we only have to find half of them.
	xl=0.5*(x2-x1);
	for (int i=0;i<m;i++) { //Loop over the desired roots.
		z=cos(3.141592653589793*(i+0.75)/(n+0.5));
		//Starting with this approximation to the ith root, we enter the main loop of refinement
		//by Newton’s method.
		do {
			p1=1.0;
			p2=0.0;
			for (int j=0;j<n;j++) { //Loop up the recurrence relation to get the
				p3=p2;				//Legendre polynomial evaluated at z.
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			//p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp; //Newton’s method.
		} while (abs(z-z1) > EPS);
		x[i]=xm-xl*z; 	  //Scale the root to the desired interval,
		x[n-1-i]=xm+xl*z; //and put in its symmetric counterpart.
		w[i]=2.0*xl/((1.0-z*z)*pp*pp); //Compute the weight
		w[n-1-i]=w[i]; //and its symmetric counterpart.
	}
}
