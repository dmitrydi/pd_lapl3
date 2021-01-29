/*
 * levin.h
 *
 *  Created on: 26 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <limits>
#include <iostream>

static const double SMALL_CONST = std::numeric_limits<double>::min()*10.0;
static const double BIG_CONST = std::numeric_limits<double>::max();
static const double EPS = std::numeric_limits<double>::epsilon();
static const double TINY = std::numeric_limits<double>::min();
static const double MAX_DOUBLE = std::numeric_limits<double>::max();
static const double LOWEST_DOUBLE = std::numeric_limits<double>::lowest();
static const int LEVIN_MAXITER = 1000;

template <typename T>
T safe(const T val) {
	if (val < LOWEST_DOUBLE) return LOWEST_DOUBLE;
	if (val > MAX_DOUBLE) return MAX_DOUBLE;
	return val;
}

struct Levin {
	//Convergence acceleration of a sequence by the Levin transformation. Initialize by calling the
	//constructor with arguments nmax, an upper bound on the number of terms to be summed, and
	//epss, the desired accuracy. Then make successive calls to the function next, which returns
	//the current estimate of the limit of the sequence. The flag cnvgd is set when convergence is
	//detected.
	std::vector<double> numer,denom; //Numerator and denominator computed via (5.3.16).
	int n,ncv;
	bool cnvgd;
	double small,big; //Numbers near machine underflow and overflow limits.
	double eps,lastval,lasteps;
	Levin(int nmax, double epss) : numer(nmax), denom(nmax), n(0), ncv(0),
			cnvgd(0), eps(epss), lastval(0.) {
		small=SMALL_CONST;
		big=BIG_CONST;
	}
	double next(double sum, double omega, double beta=1.) {
	//Arguments: sum, the nth partial sum of the sequence; omega, the nth remainder estimate
	//!n, usually from (5.3.19); and the parameter beta, which should usually be set to 1, but
	//sometimes 0.5 works better. The current estimate of the limit of the sequence is returned.
		double val;
		if (std::abs(omega) < TINY) {
			val = sum;
		} else {
			int j;
			double fact,ratio,term;
			term=1.0/(beta+n);
			denom[n]=term/omega;
			numer[n]=sum*denom[n];

			if (n > 0) {
				ratio=(beta+n-1)*term;
				for (j=1;j<=n;j++) {
					fact=(n-j+beta)*term;
					numer[n-j]=numer[n-j+1]-fact*numer[n-j];
					denom[n-j]=denom[n-j+1]-fact*denom[n-j];
					term=term*ratio;
				}
			}
			n++;
			if (std::abs(denom[0]) < small) {
				val = lastval;
			} else {
				val = numer[0]/denom[0];
			}
		}
		lasteps = std::abs(val-lastval);
		if (lasteps <= eps*std::abs(val) || lasteps <= TINY) ncv++;
		if (ncv >= 2) cnvgd = 1;
		lastval = val;
		return lastval;
	}
	//denom[188] = 4.67239e+288
	double next_verbose(double sum, double omega, double beta=1.) {
	//Arguments: sum, the nth partial sum of the sequence; omega, the nth remainder estimate
	//!n, usually from (5.3.19); and the parameter beta, which should usually be set to 1, but
	//sometimes 0.5 works better. The current estimate of the limit of the sequence is returned.
		std::cerr<<"Input: sum = " << sum << " omega = " << omega << " n = " << n << std::endl;
		double val;
		if (std::abs(omega) < TINY) {
			val = sum;
		} else {
			int j;
			double fact,ratio,term;
			term=1.0/(beta+n);
			std::cerr << "term = " << term << std::endl;
			denom[n]=term/omega;
			denom[n] = safe(denom[n]);
			numer[n]=sum*denom[n];
			numer[n] = safe(numer[n]);
			std::cerr << "denom[" << n <<"] = " << denom[n] << std::endl;
			std::cerr << "numer[" << n <<"] = " << numer[n] << std::endl;
			if (n > 0) {
				ratio=(beta+n-1)*term;
				std::cerr << "ratio = " << ratio << std::endl;
				std::cerr << "Start of cycle over j---------------" << std::endl;
				for (j=1;j<=n;j++) {
					std::cerr << "j = " << j << std::endl;
					fact=(n-j+beta)*term;
					std::cerr << "fact = " << fact << std::endl;
					numer[n-j]=numer[n-j+1]-fact*numer[n-j];
					numer[n-j] = safe(numer[n-j]);
					std::cerr << "numer[" <<n-j<<"] = " << numer[n-j] << std::endl;
					denom[n-j]=denom[n-j+1]-fact*denom[n-j];
					denom[n-j] = safe(denom[n-j]);
					std::cerr << "denom[" <<n-j<<"] = " << denom[n-j] << std::endl;
					term=term*ratio;
					std::cerr << "term = " << term << std::endl;
				}
				std::cerr << "End of cycle over j---------------" << std::endl;
			}
			n++;
			if (std::abs(denom[0]) < small) {
				val = lastval;
			} else {
				val = numer[0]/denom[0];
			}
			std::cerr << "val = " << val << std::endl;
		}
		lasteps = std::abs(val-lastval);
		std::cerr << "lastep = " << lasteps << std::endl;
		if (lasteps <= eps*std::abs(val) || lasteps <= TINY) ncv++;
		if (ncv >= 2) cnvgd = 1;
		lastval = val;
		return lastval;
	}
	double next_safe(double sum, double omega, double beta=1.) {
		//Arguments: sum, the nth partial sum of the sequence; omega, the nth remainder estimate
		//!n, usually from (5.3.19); and the parameter beta, which should usually be set to 1, but
		//sometimes 0.5 works better. The current estimate of the limit of the sequence is returned.
			double val;
			if (std::abs(omega) < TINY) {
				val = sum;
			} else {
				int j;
				double fact,ratio,term;
				term=1.0/(beta+n);
				denom[n]=safe(term/omega);
				//denom[n] = safe(denom[n]);
				numer[n]=safe(sum*denom[n]);
				//numer[n] = safe(numer[n]);
				if (n > 0) {
					ratio=(beta+n-1)*term;
					for (j=1;j<=n;j++) {
						fact=(n-j+beta)*term;
						numer[n-j]=safe(numer[n-j+1]-fact*numer[n-j]);
						//numer[n-j] = safe(numer[n-j]);
						denom[n-j]=safe(denom[n-j+1]-fact*denom[n-j]);
						//denom[n-j] = safe(denom[n-j]);
						term=term*ratio;
					}
				}
				n++;
				if (std::abs(denom[0]) < small) {
					val = lastval;
				} else {
					val = numer[0]/denom[0];
				}
			}
			lasteps = std::abs(val-lastval);
			if (lasteps <= eps*std::abs(val) || lasteps <= TINY) ncv++;
			if (ncv >= 2) cnvgd = 1;
			lastval = val;
			return lastval;
		}
};

//template <typename F>
//double LevinDSum(const F& func, const bool from_zero_i = false, const double eps = EPS, const bool verbose=false) {
//	Levin lev(LEVIN_MAXITER, eps);
//	double sum = from_zero_i ? func(0): 0.;
//	double ans;
//	const double beta = 1.;
//	int sign = 1;
//	double omega = sign*Wijn::SumNeum(func, 1);
//	double d;
//	for (int64_t r = 1; r < LEVIN_MAXITER; ++r) {
//		d = omega;
//		sum += d;
//		sign *= -1;
//		omega = sign * Wijn::SumNeum(func, r + 1);
//		ans = lev.next_safe(sum, omega, beta);
//		if (lev.cnvgd) {
//			return ans;
//		}
//
//	}
//	std::ostringstream os;
//	os << "Levin::LevinDSum did not converge in " << LEVIN_MAXITER << " steps, " << " eps = " << std::scientific << std::abs(lev.lasteps/lev.lastval);
//	throw std::runtime_error(os.str());
//}
