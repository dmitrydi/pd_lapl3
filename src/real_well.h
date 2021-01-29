/*
 * real_well.h
 *
 *  Created on: 20 ���. 2020 �.
 *      Author: Dmitry_Di
 */

#pragma once

#include "lapl_well.h"

#include <vector>

const int NCOEF = 20;

std::vector<double> CalcStehf(const int n);

class RealWell{
public:
	RealWell(const double xwd, const double xed, const double ywd, const double yed, const double Fcd, const double alpha = 0.);
	double pwd(const double td);
private:
	std::vector<double> stehf_coefs;
	LaplWell lwell;
};
