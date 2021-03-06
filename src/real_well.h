/*
 * real_well.h
 *
 *  Created on: 20 ���. 2020 �.
 *      Author: Dmitry_Di
 */

#pragma once

#include "lapl_well.h"

#include <vector>
#include <Eigen/Dense>

const int NCOEF = 12;

std::vector<double> CalcStehf(const int n);

class RealWell{
public:
	RealWell(const double xwd, const double xed, const double ywd, const double yed, const double Fcd, const double alpha = 0.);
	void PrintSourceMatrix() const;
	double pwd(const double td);
	double qwd(const double td);
	double pd(const double td, const double xd, const double yd);
	Eigen::MatrixXd pd(const double td, const Eigen::VectorXd& xd, const Eigen::VectorXd& yd);
	double pwd_(const double td);
private:
	std::vector<double> stehf_coefs;
	LaplWell lwell;
};
