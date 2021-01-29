/*
 * lapl_well.h
 *
 *  Created on: 26 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "chbessel.h"
#include "bessel.h"
#include "profile.h"

static const int NSEG = 20;
static const double SUM_EPS = 1e-14;
static const int LEVIN_MAXITER = 100;
//static const double TINY = std::numeric_limits<double>::min();
static const int KMAX = 10000;
static const int KMIN = 10;


class LaplWell {
public:
	LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_, const double alpha_ = 0.);
	double pd(const double u, const double xd, const double yd);
	double pwd(const double u);
	double qwd(const double u);

private:
	const double PI = 3.141592653589793;
	const int nseg = NSEG;
	const double eps = SUM_EPS;
	const double xwd, xed, xede, ywd, yed, Fcd, alpha;
	double iF1(const double u, const double x1, const double x2, const double yd) const; //+
	double iF2E(const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double iF2Ek(const int k, const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double i1F2H(const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double i1F2Hk(const int k, const double u, const double x1, const double x2, const double xd, const double yd, const double beta) const; //+
	double i2F2H(const double u, const double x1, const double x2,  const double yd) const;
	double SEXP(const double y, const double e) const;
	FastBessel::Bess bess;
	Eigen::MatrixXd source_matrix;
	Eigen::VectorXd rhs;
	void MakeMatrix(const double u, const double yd);
	void MakeRhs(const double u);
};
