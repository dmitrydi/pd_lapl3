/*
 * real_well.h
 *
 *  Created on: 20 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "lapl_well.h"

#include <vector>

const int NCOEF = 12;

std::vector<double> CalcStehf(const int n);

class RealWell{
public:
	RealWell(const double xwd, const double xed, const double ywd, const double yed, const double Fcd, const double alpha = 0.);
	void PrintSourceMatrix() const;
	double pwd(const double td);
	double pwd_(const double td);
private:
	std::vector<double> stehf_coefs;
	LaplWell lwell;
};
