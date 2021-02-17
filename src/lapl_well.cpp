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
		xwd(xwd_), xed(xed_), xede(xed_), ywd(ywd_), yed(yed_), Fcd(Fcd_), alpha(alpha_),dx(1./NSEG), bess(false),
		source_matrix(2*NSEG+1, 2*NSEG+1),
		rhs(2*NSEG+1), _src_matrix(2*NSEG, 2*NSEG),
		if1_matrix(2*NSEG, 2*NSEG), if2e_matrix(2*NSEG, 2*NSEG), i1f2h_matrix(2*NSEG, 2*NSEG),
		i2f2h_matrix(2*NSEG,2*NSEG),  i1f2h_buf(2*NSEG), if2e_buf(2*NSEG){
	_src_matrix = MakeSrcMatrix();
};

double LaplWell::pwd_(const double u) {
	MakeMatrix(u, ywd);
	MakeRhs(u);
	return source_matrix.colPivHouseholderQr().solve(rhs)(0);
}

double LaplWell::qwd(const double u) {
	return pwd(u)/u;
}

double LaplWell::pwd(const double u) {
	NewMakeMatrix(u);
	MakeRhs(u);
	return source_matrix.colPivHouseholderQr().solve(rhs)(0);
}

double LaplWell::pd(const double u, const double xd, const double yd) {
	NewMakeMatrix(u);
	MakeRhs(u);
	auto svect = source_matrix.colPivHouseholderQr().solve(rhs);
	auto green = MakeGreenVector(u, xd, yd);
	double ans = 0.;
	for (int i = 0; i < 2*NSEG; ++i) {
		ans += svect(i+1)*green(i);
	}
	return ans;
}

Eigen::MatrixXd LaplWell::pd(const double u, const Eigen::VectorXd& xd, const Eigen::VectorXd& yd) {
	size_t nxd, nyd;
	nxd = xd.size();
	nyd = yd.size();
	Eigen::MatrixXd ans;
	ans = Eigen::MatrixXd::Zero(nxd, nyd);
	NewMakeMatrix(u);
	MakeRhs(u);
	Eigen::VectorXd svect = source_matrix.colPivHouseholderQr().solve(rhs);
	for (int i = 0; i < nxd; ++i) {
		for (int j = 0; j< nyd; ++j) {
			auto green = MakeGreenVector(u, xd(i), yd(j));
			for (int k = 0; k < 2*NSEG; ++k) {
				ans(i, j) += svect(k+1)*green(k);
			}
		}
	}
	return ans;
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
	double d, dd, ek;
	double A;
	for (int k = 1; k < KMAX; ++k) {
		dd = iF2Ek(k, u, x1, x2, xd, yd);
		sum += dd;
		ek = sqrt(u+k*k*PI*PI/xed/xed+alpha*alpha);
		A = 2*xed/PI/(1-exp(-2*ek*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (isnan(d)) d = 0.;
		if (k > KMIN && (abs(d) <= TINY || abs(sum) <= TINY || abs(d/sum) < eps)) {
			return sum;
		}
	}
	double eps = abs(d/sum);
	ostringstream os;
	os << "iF2E did not converge in 10000 steps, last eps = " << scientific << eps  << " sum = " << sum<< endl;
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
			if (isnan(d)) d = 0.;
			if (k>kmin && (d <= TINY || sum <= TINY || abs(d/sum) < eps)) {
				return 0.5*xede/PI*sum;

			}
		}
		double eps = abs(d/sum);
		ostringstream os;
		os << "i1F2H did not converge in 10000 steps, last eps = " << scientific << eps  << " sum = " << sum<< endl;
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
			return mult*(bess.ik00x(squ*abs(t1))+bess.ik00x(squ*abs(t2)));
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

Eigen::MatrixXd LaplWell::MakeSrcMatrix() const {
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			ans(i,j) = 0.;
		}
	}
	double dx = 1./NSEG;
	double dx2_8 = 0.125*dx*dx;
	double dx2_2 = 0.5*dx*dx;
	double coef = PI/Fcd;
	for (int j = 0; j < NSEG; ++j) {
		ans(j+NSEG,j+NSEG) = coef*dx2_8;
		ans(NSEG-j-1, NSEG-j-1) = ans(j+NSEG,j+NSEG);
		double xj = dx*(j+0.5);
		for (int i = 0; i < j; ++i) {
			ans(j+NSEG,i+NSEG) = coef*(dx2_2 + dx*(xj - (i+1)*dx));
			ans(NSEG-j-1, NSEG-i-1) = ans(j+NSEG,i+NSEG);
		}
	}
	return ans;
}

void LaplWell::PrintSourceMatrix() const {
	auto ans = MakeSrcMatrix();
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			cout << fixed << setprecision(7) << ans(i, j) << ' ';
		}
		cout << endl;
	}
}



void LaplWell::MakeMatrix(const double u, const double yd) {
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	for (int i = 0; i < 2*nseg; ++i) {
		source_matrix(i,0) = 1.;
		source_matrix(2*nseg, i+1) = 1.;
	}
	//auto src_m = MakeSrcMatrix();
	for (int i = 0; i< 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			source_matrix(i, j+1) = _src_matrix(i,j);
		}
	}
	// double iF1(const double u, const double x1, const double x2, const double yd) const
	{ //LOG_DURATION("iF1");
	for (int i = 0; i < 2*nseg; ++i) {
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*iF1(u, x1, x2, yd);
		}
	}
	}
	// double iF2E(const double u, const double x1, const double x2, const double xd, const double yd) const
	{//LOG_DURATION("iF2E");
	for (int i = 0; i < 2*nseg; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += 1.*mult*iF2E(u, x1, x2, xd, yd); //?
		}
	}
	}
	// double i1F2H(const double u, const double x1, const double x2, const double xd, const double yd) const
	{//LOG_DURATION("i1F2H");
	for (int i = 0; i < 2*nseg; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*i1F2H(u, x1, x2, xd, yd);
		}
	}
	}
	// double i2F2H(const double u, const double x1, const double x2,  const double yd) const
	{//LOG_DURATION("i2F2H");
	for (int i = 0; i < 2*nseg; ++i) {
		for (int j = 0; j < 2*nseg; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			source_matrix(i, j+1) += mult*i2F2H(u, x1, x2, yd);
		}
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


// ----new mtrix versions
//
//
// ----------------------

void LaplWell::NewMakeMatrix(const double u) {
	double mult = -1.*PI/xed;
	for (int i = 0; i < 2*nseg; ++i) {
		source_matrix(i,0) = 1.;
		source_matrix(2*nseg, i+1) = 1.;
	}
	fill_if1(u);
	fill_if2e(u);
	fill_i1f2h(u);
	fill_i2f2h(u);
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			source_matrix(i, j+1) = mult*(if1_matrix(i,j)+ if2e_matrix(i,j)+ i1f2h_matrix(i,j)+ i2f2h_matrix(i,j)) + _src_matrix(i,j);
		}
	}
}

Eigen::VectorXd LaplWell::MakeGreenVector(const double u, const double xd, const double yd) {
	Eigen::VectorXd ans(2*NSEG), buf(2*NSEG);
	double mult = PI/xed;
	vect_if1_yd(u, yd, buf);
	ans = mult*buf;
	vect_if2e_yd(u, xd, yd, buf);
	ans += mult*buf;
	vect_i1f2h(u, xd, buf);
	ans += mult*buf;
	vect_i2f2h_yd(u, yd, buf);
	ans += mult*buf;
	return ans;
}

Eigen::MatrixXd LaplWell::show_matrix(const double u) {
	return show_if1_matrix(u)+show_i1f2h_matrix(u)+show_i2f2h_matrix(u)+show_if2e_matrix(u);
}

void LaplWell::fill_i2f2h(const double u) {
	double squ = sqrt(u+alpha*alpha);
	i2f2h_matrix = -0.5*exp(-squ*abs(ywd-ywd))/squ*(dx)*Eigen::MatrixXd::Ones(2*NSEG, 2*NSEG);
}

void LaplWell::vect_i2f2h_yd(const double u, const double yd, Eigen::VectorXd& buf) {
	double squ = sqrt(u+alpha*alpha);
	buf = -0.5*exp(-squ*abs(yd-ywd))/squ*(dx)*Eigen::VectorXd::Ones(2*NSEG);
}

void LaplWell::fill_if1(const double u) {
	double squ = sqrt(u+alpha*alpha);
	double ans = 0.5*dx/squ;
	double dy = std::abs(ywd-ywd);
	double sumy = ywd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+SEXP(yed, squ));
	if1_matrix = Eigen::MatrixXd::Ones(2*NSEG, 2*NSEG)*ans;
}

void LaplWell::vect_if1_yd(const double u, const double yd, Eigen::VectorXd& buf) {
	double squ = sqrt(u+alpha*alpha);
	double ans = 0.5*dx/squ;
	double dy = std::abs(yd-ywd);
	double sumy = yd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+SEXP(yed, squ));
	buf = Eigen::VectorXd::Ones(2*NSEG)*ans;
}

void LaplWell::fill_if2e(const double u) {
	if2e_matrix = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	static const double term1 = PI/xed*(2*yed-(ywd+ywd));
	static const double term2 = PI/xed*(2*yed-abs(ywd-ywd));
	static const double term3 = PI/xed*(ywd+ywd);
	static const double term4 = PI/xed*(2*yed+abs(ywd-ywd));
	static const double dterm1 = 1.-exp(-term1);
	static const double dterm2 = 1.-exp(-term2);
	static const double dterm3 = 1.-exp(-term3);
	static const double dterm4 = 1.-exp(-term4);
	double ek_term, ek_, sexp_, aydywd, kpiOxed, ydPywd, mmult, A, d, max_mat=0.;
	for (int k = 1; k <= KMAX; ++k) {
		ek_term = k*PI/xede;
		ek_ = sqrt(u + ek_term*ek_term + alpha*alpha);
		sexp_ = SEXP(yed, ek_);
		aydywd = abs(ywd-ywd); //!
		kpiOxed = k*PI/xed;
		ydPywd = ywd + ywd; //!
		mmult = 2./kpiOxed/ek_*((std::exp(-ek_*(2.*yed - ydPywd)) + std::exp(-ek_*ydPywd) + std::exp(-ek_*(2.*yed-aydywd)))*(1 + sexp_)
				+ std::exp(-ek_*aydywd)*sexp_);
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1.+j*dx;
			double x2 = x1 + dx;
			if2e_buf(j)  = mmult*sin(0.5*kpiOxed*(x2 - x1))*std::cos(0.5*kpiOxed*(2.*xwd + x1 + x2));
			if (abs(if2e_buf(j)) > max_mat) max_mat = abs(if2e_buf(j));
		}
		for (int i = 0; i < 2*NSEG; ++i) {
			double xd = (xwd-1.+(i+0.5)*dx);
			double row_mult = cos(kpiOxed*xd);
			for (int j = 0; j < 2*NSEG; ++j) {
				if2e_matrix(i,j) += row_mult*if2e_buf(j);
			}
		}
		A = 2*xed/PI/(1-exp(-2*ek_*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (isnan(d)) d = 0.;
		if (k > KMIN && (abs(d) <= TINY || abs(max_mat) <= TINY || abs(d/max_mat) < eps)) break;
	}
}


void LaplWell::vect_if2e_yd(const double u, const double xd, const double yd, Eigen::VectorXd& buf) {
	buf = Eigen::VectorXd::Zero(2*NSEG);
	static const double term1 = PI/xed*(2*yed-(yd+ywd));
	static const double term2 = PI/xed*(2*yed-abs(yd-ywd));
	static const double term3 = PI/xed*(yd+ywd);
	static const double term4 = PI/xed*(2*yed+abs(yd-ywd));
	static const double dterm1 = 1.-exp(-term1);
	static const double dterm2 = 1.-exp(-term2);
	static const double dterm3 = 1.-exp(-term3);
	static const double dterm4 = 1.-exp(-term4);
	double ek_term, ek_, sexp_, aydywd, kpiOxed,ydPywd, mmult, A, d, max_mat=0., row_mult;
	for (int k = 1; k <= KMAX; ++k) {
		ek_term = k*PI/xede;
		ek_ = sqrt(u + ek_term*ek_term + alpha*alpha);
		sexp_ = SEXP(yed, ek_);
		aydywd = abs(yd-ywd); //!
		kpiOxed = k*PI/xed;
		ydPywd = yd + ywd; //!
		mmult =  2./kpiOxed/ek_*((std::exp(-ek_*(2.*yed - ydPywd)) + std::exp(-ek_*ydPywd) + std::exp(-ek_*(2.*yed-aydywd)))*(1 + sexp_)
				+ std::exp(-ek_*aydywd)*sexp_);
		row_mult = cos(kpiOxed*xd);
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1.+j*dx;
			double x2 = x1 + dx;
			buf(j)  += row_mult*mmult*sin(0.5*kpiOxed*(x2 - x1))*std::cos(0.5*kpiOxed*(2.*xwd + x1 + x2));
			if (abs(buf(j)) > max_mat) max_mat = abs(buf(j));
		}
		A = 2*xed/PI/(1-exp(-2*ek_*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (isnan(d)) d = 0.;
		if (k > KMIN && (abs(d) <= TINY || abs(max_mat) <= TINY || abs(d/max_mat) < eps)) break;
	}
}

void LaplWell::fill_i1f2h(const double u) {
	for (int i: {0, 2*NSEG-1}) {
		double xd = xwd-1.+(i+0.5)*dx;
		vect_i1f2h(u, xd, i1f2h_buf);
		for (int j = 0; j < 2*NSEG; ++j)
			i1f2h_matrix(i,j) = i1f2h_buf(j);
	}
	for (int i = 1; i < 2*NSEG - 1; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			if (j >= i) {
				i1f2h_matrix(i,j) = i1f2h_matrix(0, j-i);
			}
			else {
				i1f2h_matrix(i,j) = i1f2h_matrix(2*NSEG-1, 2*NSEG-1-i+j);
			}
		}
	}
}

void LaplWell::vect_i1f2h(const double u, const double xd, Eigen::VectorXd& buf) {
	buf = Eigen::VectorXd::Zero(2*NSEG);
	const double squ = sqrt(u+alpha*alpha);
	static const int kmin = static_cast<int>(abs(0.5*(2./squ/xede-(-1.+(2*NSEG-1)*dx)/xed+(xwd-1+0.5*dx)/xed+xwd/xed)));
	static const double dmult = 8.*dx*exp(squ*xede/xed*abs((xwd-1+0.5*dx) + xwd - (-1.+(2*NSEG-1)*dx)))/1.-exp(-squ*2.*xede); //?
	double mult = xed/xede/squ*0.5*xede/PI;
	double t1, t2, x1, x2, elem, d;
	for (int k = 0; k <=KMAX; ++k) {
		for (double beta: {-1., 1.}) {
			for (int j = 0; j < 2*NSEG; ++j) {
				x1 = -1.+j*dx;
				x2 = x1 + dx;
				t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed-2.*k);
				t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed-2.*k);
				elem = mult*bess.abs_ik0ab(t1, t2);
				if (k > 0) {
					t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed+2.*k);
					t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed+2.*k);
					elem += mult*bess.abs_ik0ab(t1, t2);
				}
				buf(j) += elem;
			}
		}
		d = dmult*exp(-squ*2*k*xede);
		if (isnan(d)) d = 0.;
		if (k>kmin && (d <= TINY || buf(2*NSEG-1) <= TINY || abs(d/buf(2*NSEG-1)) < eps)) break;
	}
}

void LaplWell::vect_i1f2h_yd(const double u, const double xd, const double yd, Eigen::VectorXd& buf) {
	double adyd = abs(ywd-yd);
	if (adyd < 1e-16)
	{
		vect_i1f2h(u, xd, buf);
	} else {
		buf = Eigen::VectorXd::Zero(2*NSEG);
		const double squ = sqrt(u+alpha*alpha);
		auto func = [this, adyd, squ](double x) {return this->bess.k0(squ*std::sqrt(x*x+adyd*adyd));};
		double t1, t2, x1, x2, elem;
		double mult = xed/xede;
		for (int k = 0; k <=KMAX; ++k) {
			for (double beta: {-1., 1.}) {
				for (int j = 0; j < 2*NSEG; ++j) {
					x1 = -1.+j*dx;
					x2 = x1 + dx;
					t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed-2.*k);
					t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed-2.*k);
					elem = mult*qromb(func, t1, t2, eps);
					if (k > 0) {
						t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed+2.*k);
						t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed+2.*k);
						elem += mult*qromb(func, t1, t2, eps);
					}
					buf(j) += elem;
				}
			}
			if ( elem <= TINY || buf(2*NSEG-1) <= TINY || abs(elem/buf(2*NSEG-1)) < eps) break;
		}
	}
}

// ---- new show
Eigen::VectorXd LaplWell::show_green_vector(const double u, const double xd, const double yd) {
	return MakeGreenVector(u, xd, yd);
}

Eigen::MatrixXd LaplWell::show_if1_matrix(const double u) {//fast
	fill_if1(u);
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	double mult = -1.*PI/xed;
	return if1_matrix*mult;
}


Eigen::MatrixXd LaplWell::show_i1f2h_matrix(const double u) {
	fill_i1f2h(u);
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	double mult = -1.*PI/xed;
	for (int i =0; i < 2*NSEG; ++i)
		for (int j=0; j<2*NSEG; ++j)
			ans(i,j) = mult*i1f2h_matrix(i,j);
	return ans;
}

Eigen::MatrixXd LaplWell::show_i2f2h_matrix(const double u) {
	fill_i2f2h(u);
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	double mult = -1.*PI/xed;
	for (int i =0; i < 2*NSEG; ++i)
		for (int j=0; j<2*NSEG; ++j)
			ans(i,j) = mult*i2f2h_matrix(i,j);
	return ans;
}

Eigen::MatrixXd LaplWell::show_if2e_matrix(const double u) {//fast
	fill_if2e(u);
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	double mult = -1.*PI/xed;
	return if2e_matrix*mult;
}



// --- old show

Eigen::MatrixXd LaplWell::make_if1(const double u) {
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	ans = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		//double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			ans(i, j) += mult*iF1( u,  x1,  x2, ywd);
		}
	}
	return ans;
}

Eigen::MatrixXd LaplWell::make_i1f2h(const double u) {
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	ans = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			ans(i, j) += mult*i1F2H(u, x1, x2, xd, ywd);
		}
	}
	return ans;
}

Eigen::MatrixXd LaplWell::make_i2f2h(const double u) {
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	ans = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			ans(i, j) += mult*i2F2H(u, x1, x2, ywd);
		}
	}
	return ans;
}

Eigen::MatrixXd LaplWell::make_if2e(const double u){
	double dx = 1./nseg;
	double mult = -1.*PI/xed;
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	ans = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		double xd = xwd-1. + (i + 0.5)*dx;
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1. + j*dx;
			double x2 = x1 + dx;
			ans(i, j) += mult*iF2E(u, x1, x2, xd, ywd);
		}
	}
	return ans;
};
