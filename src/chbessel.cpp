/*
 * chbessel.cpp
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */
#include "chbessel.h"

using namespace std;

namespace FastBessel {

vector<long double> Bess::_coef() {
	vector<long double> ans(MAXIT_IKBESS);
	long double dum;
	long double inv_k_fact = 1.;
	for (int i = 1; i < MAXIT_IKBESS; ++i) {
		dum = 1./(long double)i;
		inv_k_fact *= dum;
		ans[i] = 1./(2.*i + 1.)*inv_k_fact*inv_k_fact;
	}
	return ans;
}

vector<long double> Bess::_ns() {
	vector<long double> ans(MAXIT_IKBESS);
	long double inv_n_sum = 0.;
	long double dum;
	for (int i = 1; i < MAXIT_IKBESS; ++i) {
		dum = 1./(long double)i;
		inv_n_sum += dum;
		ans[i] = inv_n_sum;
	}
	return ans;
}

Bess::Bess(const int n, const int m, const double dd): d(dd), N(n), M(m),
		ak(fak(N, N, d)), ck(fck(M, M, d, 0.)), cha(N+1, 1.), chc(M+1, 1.),
		coef(_coef()), ns(_ns()), bess(), gs(GAUSS_POINTS) {};

double Bess::num_ik0ab(const double x1, const double x2) const {
	return gs.Integrate(bess, x1, x2);
}


long double Bess::_ik02(const double x) const {
	if (x < 0.) throw;
	if (x <= TINY) return 0.;
	long double xO2=0.5*x;
	const long double A = -(log(xO2) + EUL_GAMMA_LD);
	long double xO2sq = xO2*xO2;
	long double sum = 0.;
	for (int i = MAXIT_IKBESS-1; i>= 1; i--) {
		sum = xO2sq*(coef[i]*(A+1./(2.*i+1.)+ns[i]) + sum);
	}
	return x*A+x+x*sum;
}

double Bess::ik02(const double x) const {
	if (x == 0.) return PI2;
	return PI2 - _ik02(x);
}

double Bess::_k0(const double x) {
	double z = d/x;
	double z2 = 2.*(2.*z-1.);
	cha[1] = (2.*z-1.);
	for (int i = 2; i<=N; ++i) {
		cha[i] = z2*cha[i-1] - cha[i-2];
	}
	double sum = 0.;
	for (int i = N; i >= 0; --i) {
		sum += cha[i]*ak[i];
	}
	return sum;
}

//double Bess::_ik0(const double x) {
//	if (x < d) throw;
//	double z = d/x;
//	double z2 = 2.*(2*z-1.);
//	chc[0] = 1;
//	chc[1] = (2.*z-1.);
//	for (int i = 2; i <= M; ++i) {
//		chc[i] = z2*chc[i-1] - chc[i-2];
//	}
//	double sum = 0.;
//	for (int i = M; i>=0; --i) {
//		sum += chc[i]*ck[i];
//	}
//	return sum;
//}

double Bess::_ik0(const double x) const {
	if (x < d) throw;
	double z = d/x;
	double z2 = 2.*(2*z-1.);
	double ch_0 = 1.;
	double ch_1 = (2.*z-1.);
	double ch_2;
	double sum = ch_0*ck[0] + ch_1*ck[1];
	for (int i = 2; i <=M; ++i) {
		ch_2 = z2*ch_1 - ch_0;
		sum += ch_2*ck[i];
		ch_0 = ch_1;
		ch_1 = ch_2;
	}
	return sum;
}

double Bess::k0(const double x) {
	return exp(-x)/sqrt(x)*_k0(x);
}

double Bess::ik0(const double x) const {
	if (x < 2.) return ik02(x);
	return exp(-x)/sqrt(x)*_ik0(x);
}

double Bess::ik0ab(const double x1, const double x2) const {
	if (x1 <= 0.1 || (x2-x1) > 0.5) {
		return analytic_ik0ab(x1, x2);
	}
	return num_ik0ab(x1, x2);
}

double Bess::analytic_ik0ab(const double x1, const double x2) const {
	if (x2 < 2.) {
		return _ik02(x2) - _ik02(x1);
	} else if (x1 < 2. && x2 >=2.) {
		return  (PI2 - _ik02(x1)) - ik0(x2) ;
	} else {
		return ik0(x1) - ik0(x2);
	}
}

vector<double> Bess::fak(const int m, const int n, const double d) {
	vector<__float128> fq = fakq(m,n,(__float128)d, SQRT_PI2Q); // here scaling is used SQRT_PI2Q
	vector<double> ans(fq.size());
	for (size_t i = 0; i< fq.size(); ++i) {
		ans[i] = (double)fq[i];
	}
	return ans;
}

vector<double> Bess::fck(const int m, const int n, const double d, const double mu) {
	vector<__float128> fq = fckq(m,n,(__float128)d, (__float128)mu, SQRT_PI2Q); // here scaling is used
	vector<double> ans(fq.size());
		for (size_t i = 0; i< fq.size(); ++i) {
			ans[i] = (double)fq[i];
		}
	return ans;
}

vector<__float128> Bess::fakq(const int m, const int n, const __float128 d, const __float128 multt) {
	vector<__float128> cfa(n+4, 0.q);
	vector<__float128> ans(m+1);
	cfa[n] = 1.q;
	__float128 cf1, cf2, cf3;
	__float128 ek;
	for (int k = n-1; k>=0; k--) {
		cf1 = 2.q*(k+1.q)*(1.q-(2.q*k+3.q)*(k+1.5q)*(k+1.5q)/(2.q*(k+2.q)*(k+0.5q)*(k+0.5q))-4.q*d/((k+0.5q)*(k+0.5q)));
		cf2 = (1.q-2.q*(k+1.q)*(2.q*k+3.q-4.q*d)/((k+0.5q)*(k+0.5q)));
		cf3 = -(k+1.q)*(k+2.5q)*(k+2.5q)/((k+2.q)*(k+0.5q)*(k+0.5q));
		ek = k==0 ? 1.q: 2.q;
		cfa[k] = ek/2.q*(cfa[k+1]*cf1+cfa[k+2]*cf2+cfa[k+3]*cf3);
	}
	__float128 un=0.q;
	__float128 mult = n%2==0? 1.0q : -1.0q;
	for (int k=n; k>=0; k--) {
		un += mult*cfa[k];
		mult *= -1.q;
	}
	for (int k=0; k<=m; k++) {
		ans[k] = multt*cfa[k]/un;
	}
	return ans;
}

vector<__float128> Bess::fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt) {
	vector<__float128> ak = fakq(n+3, n+3, d);
	vector<__float128> dk(n+4, 0.q);
	vector<__float128> eek(n+4, 0.q);
	__float128 ek = 2.q;
	dk[n] = -ek/2.q/(n+0.5q-mu)*4.q*d*(ak[n+2]-ak[n+1]); // ?
	for (int k=n-1; k>=0; k--) {
		if (k==0) ek=1.q;
		dk[k] = -ek/2.q/(k+0.5q-mu)*(4.q*d*(ak[k+2]-ak[k+1]) +
				(3.q*k+3.5q-mu+4.q*d)*dk[k+1] +
				(3.q*k+5.5q+mu-4.q*d)*dk[k+2] +
				(k+2.5q+mu)*dk[k+3]);
		if (k==n-1) dk[n] = 0.; // ?
	}
	ek = 2.0q;
	eek[n] = 1.q;
	for (int k=n-1; k>=0; k--) {
		if (k==0) ek=1.q;
		eek[k] = -ek/2.q/(k+0.5q-mu)*( (3.q*k+3.5q-mu+4.q*d)*eek[k+1] + (3.q*k+5.5q+mu-4.q*d)*eek[k+2] +
				(k+2.5q+mu)*eek[k+3]);
	}
	__float128 sumd=0.q, sume=0.q;
	__float128 mult = n%2==0? 1.0q : -1.0q;
	for (int k=n; k>=0; k--) {
		sume += mult*eek[k];
		sumd += mult*dk[k];
		mult *= -1.q;
	}
	__float128 An = (1.q-sumd)/sume;
	vector<__float128> ans(m+1);
	for (int k=0; k<=m; ++k) {
		ans[k] = multt*(dk[k]+An*eek[k]);
	}
	return ans;
}

}


