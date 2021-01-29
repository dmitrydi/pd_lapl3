/*
 * lapl_well.cpp
 *
 *  Created on: 26 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#include "lapl_well.h"

using namespace std;
using namespace FastBessel;

LaplWell::LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_, const double alpha_):
		xwd(xwd_), xed(xed_), xede(xed_), ywd(ywd_), yed(yed_), Fcd(Fcd_), alpha(alpha_), source_matrix(2*NSEG+1, 2*NSEG+1),
		rhs(2*NSEG+1) {};

double LaplWell::pwd(const double u) {
	MakeMatrix(u, ywd);
	MakeRhs(u);
	return source_matrix.colPivHouseholderQr().solve(rhs)(0);
}

double LaplWell::iF1(const double u, const double x1, const double x2, const double yd) const {
	double squ = sqrt(u+alpha*alpha);
	double ans = 0.5*(x2-x1)/squ;
	double dy = std::abs(yd-ywd);
	double sumy = yd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+SEXP(yed, squ));
	return ans;
}

double LaplWell::iF2E(const double u, const double x1, const double x2, const double xd, const double yd) const {
	double term1 = PI/xed*(2*yed-(yd+ywd));
	double term2 = PI/xed*(2*yed-abs(yd-ywd));
	double term3 = PI/xed*(yd+ywd);
	double term4 = PI/xed*(2*yed+abs(yd-ywd));
	double dterm1 = 1.-exp(-term1);
	double dterm2 = 1.-exp(-term2);
	double dterm3 = 1.-exp(-term3);
	double dterm4 = 1.-exp(-term4);
	double sum = 0.;
	double d, ek;
	double A;
	for (int k = 1; k < KMAX; ++k) {
		sum += iF2Ek(k, u, x1, x2, xd, yd);
		ek = sqrt(u+k*k*PI*PI/xed/xed+alpha*alpha);
		A = 2*xed/PI/(1-exp(-2*ek*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (k > KMIN && (abs(d/sum) < eps)) {
			return sum;
		}
	}
	double eps = abs(d/sum);
	ostringstream os;
	os << "iF2e did not converge in 10000 steps, last eps = " << scientific << eps  << " " << sum<< endl;
	throw runtime_error(os.str());
}

double LaplWell::iF2Ek(const int k, const double u, const double x1, const double x2, const double xd, const double yd) const {
	double ek_term = k*PI/xede;
	double ek_ = sqrt(u + ek_term*ek_term + alpha*alpha);
	double sexp_ = SEXP(yed, ek_);
	double aydywd = std::abs(yd-ywd);
	double kpiOxed = k*PI/xed;
	double ydPywd = yd + ywd;
	double res;
	res = 2./kpiOxed/ek_;
	res *= std::cos(kpiOxed*xd);
	res *= std::sin(0.5*kpiOxed*(x2 - x1));
	res *= std::cos(0.5*kpiOxed*(2.*xwd + x1 + x2));
	res *= (std::exp(-ek_*(2.*yed - ydPywd)) + std::exp(-ek_*ydPywd) + std::exp(-ek_*(2.*yed-aydywd)))*(1 + sexp_)
			+ std::exp(-ek_*aydywd)*sexp_;
	return res;
}

double LaplWell::i1F2H(const double u, const double x1, const double x2, const double xd, const double yd) const {
	// assure sqrt(u+a^2)*abs(t) >= 2
		double beta1 = 1., beta2 = -1.;
		double squ = sqrt(u+alpha*alpha);
		int kmin = static_cast<int>(abs(0.5*(2./squ/xede-x1/xed+xd/xed+xwd/xed)));
		kmin = kmin > 1 ? kmin : 1;
		double sum = i1F2Hk(0, u, x1, x2, xd, yd, beta1)+i1F2Hk(0, u, x1, x2, xd, yd, beta2);
		double d;
		double dd;
		double dlm = 1.-exp(-squ*2.*xede);
		double dt = abs(x1-x2)*xede/xed;
		double mult = 8.*xed*dt/xede*exp(squ*xede*abs(xd/xed+xwd/xed-x1/xed));
		for (int k = 1; k <=KMAX; ++k) {
			dd= i1F2Hk(k, u, x1, x2, xd, yd, beta1);
			dd+= i1F2Hk(-k, u, x1, x2, xd, yd, beta1);
			dd += i1F2Hk(k, u, x1, x2, xd, yd, beta2);
			dd += i1F2Hk(-k, u, x1, x2, xd, yd, beta2);
			sum += dd;
			d = mult*exp(-squ*2*k*xede)/dlm;
			if (k>kmin && (abs(d/sum) < eps || abs(sum) < TINY)) {
				return 0.5*xede/PI*sum;
			}
		}
		double eps = abs(d/sum);
		ostringstream os;
		os << "i1F2H did not converge in 10000 steps, last eps = " << scientific << eps  << " " << sum<< endl;
		throw runtime_error(os.str());
}

double LaplWell::i1F2Hk(const int k, const double u, const double x1, const double x2, const double xd, const double yd, const double beta) const {
	double t = xede*(xd/xed+beta*xwd/xed-2*k);
	double t1 = t - xede*x2/xed;
	double t2 = t - xede*x1/xed;
	double adyd = std::abs(yd-ywd);
	if (adyd <= TINY) {
		double squ = sqrt(u+alpha*alpha);
		double mult = xed/xede/squ;
		if (t1 >= 0.0) {
			return mult*bess.ik0ab(squ*t1, squ*t2);
		} else if (t2 <= 0.) {
			return mult*bess.ik0ab(squ*abs(t2), squ*abs(t1));
		} else {
			return mult*(bess.ik0ab(0., squ*abs(t1))+bess.ik0ab(0., squ*abs(t2)));
		}
	} else {
		throw logic_error("i1F2Hk for yd != ywd not implemented");
	}

}

double LaplWell::i2F2H(const double u, const double x1, const double x2,  const double yd) const {
	double squ = sqrt(u+alpha*alpha);
	return -0.5*exp(-squ*abs(yd-ywd))/squ*(x2-x1);
}

double LaplWell::SEXP(const double y, const double e) const {
	double b = exp(-2.*y*e);
	return b/(1.-b);
}

void LaplWell::MakeMatrix(const double u, const double yd) {
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	for (int i = 0; i < 2*nseg; ++i) {
		source_matrix(i,0) = 1.;
		source_matrix(2*nseg, i+1) = 1.;
	}
	// double iF1(const double u, const double x1, const double x2, const double yd) const
	for (int i = 0; i < 2*nseg; ++i) {
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) = mult*iF1(u, x1, x2, yd);
		}
	}
	// double iF2E(const double u, const double x1, const double x2, const double xd, const double yd) const
	for (int i = 0; i < 2*nseg; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*iF2E(u, x1, x2, xd, yd);
		}
	}
	// double i1F2H(const double u, const double x1, const double x2, const double xd, const double yd) const
	for (int i = 0; i < 2*nseg; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*i1F2H(u, x1, x2, xd, yd);
		}
	}
	// double i2F2H(const double u, const double x1, const double x2,  const double yd) const
	for (int i = 0; i < 2*nseg; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*i2F2H(u, x1, x2, yd);
		}
	}
}

void LaplWell::MakeRhs(const double u) {
	double coef = PI/Fcd/nseg/u;
	for (int i = 0; i < nseg; ++i) {
		rhs(nseg+i) = coef*(i+0.5);
		rhs(nseg-i-1) = rhs[nseg+i];
	}
	rhs(2*nseg) = 2.*nseg/u;
}

