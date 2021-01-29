/*
 * real_well.cpp
 *
 *  Created on: 20 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "real_well.h"

std::vector<double> CalcStehf(const int n) { //OK
	std::vector<double> v(n+1);
	std::vector<double> g(161);
	std::vector<double> h(81);
	g[1] = 1.;
	int NH = n/2;
	for (int i = 2; i < n+1; ++i) {
		g[i] = g[i-1]*i;
	}
	h[1] = 2./g[NH-1];
	for (int i = 2; i < NH+1; i++){
		double fi = i;
		if (i != NH) {
			h[i] = std::pow(fi, NH)*g[2*i]/(g[NH - i]*g[i]*g[i - 1]);
		} else {
			h[i] = std::pow(fi, NH)*g[2*i]/(g[i]*g[i - 1]);
		}
	}
	int SN = 2 * (NH - (NH/2)*2) - 1;
	for (int i = 1; i < n+1; ++i) {
		v[i] = 0.;
		int K1 = (i+1)/2;
		int K2 = i;
		if(K2>NH) K2=NH;
		for (int k = K1; k < K2+1; ++k) {
			if (2*k - i == 0) {
				v[i] += h[k]/g[i-k];
			} else if (i==k) {
				v[i] += h[k]/g[2*k - i];
			} else {
				v[i] += h[k]/(g[i - k]*g[2*k - i]);
			}
		}
		v[i] *= SN;
		SN = -1*SN;
	}
	return v;
}

RealWell::RealWell(const double xwd, const double xed, const double ywd, const double yed, const double Fcd, const double alpha):
		stehf_coefs(CalcStehf(NCOEF)),
		lwell(LaplWell(xwd,xed,ywd,yed,Fcd, alpha)){
};
double RealWell::pwd(const double td) {
	double s_mult = std::log(2.)/td;
	double ans = 0.;
	for (int i = 1; i <=NCOEF; ++i) {
		double s = i*s_mult;
		ans += lwell.pwd(s)*s*stehf_coefs[i]/i;
	}
	return ans;
};
