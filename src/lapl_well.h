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
#include "profile.h"
#include "quadrature.h"

static const int NSEG = 20;
static const double SUM_EPS = 1e-8;
static const int LEVIN_MAXITER = 100;
//static const double TINY = std::numeric_limits<double>::min();
static const int KMAX = 10000;
static const int KMIN = 10;


class LaplWell {
public:
	LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_, const double alpha_ = 0.);
	void PrintSourceMatrix() const;
	double pd(const double u, const double xd, const double yd);
	Eigen::MatrixXd pd(const double u, const Eigen::VectorXd& xd, const Eigen::VectorXd& yd);
	double pwd(const double u);
	double pwd_(const double u);
	double qwd(const double u);
	Eigen::MatrixXd make_i1f2h(const double u); // slow
	Eigen::MatrixXd make_i2f2h(const double u); // slow
	Eigen::MatrixXd make_if2e(const double u); // slow
	Eigen::MatrixXd make_if1(const double u);
	Eigen::MatrixXd show_i1f2h_matrix(const double u); //fast OK
	Eigen::MatrixXd show_i2f2h_matrix(const double u); //fast OK
	Eigen::MatrixXd show_if2e_matrix(const double u); //fast OK
	Eigen::MatrixXd show_if1_matrix(const double u);
	Eigen::MatrixXd show_matrix(const double u);
	Eigen::VectorXd show_green_vector(const double u, const double xd, const double yd);

private:
	const double PI = 3.141592653589793;
	const int nseg = NSEG;
	const double eps = SUM_EPS;
	const double xwd, xed, xede, ywd, yed, Fcd, alpha, dx;
	FastBessel::Bess bess;
	double iF1(const double u, const double x1, const double x2, const double yd) const; //+
	double iF2E(const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double iF2Ek(const int k, const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double i1F2H(const double u, const double x1, const double x2, const double xd, const double yd) const; //+
	double i1F2Hk(const int k, const double u, const double x1, const double x2, const double xd, const double yd, const double beta) const; //+
	double i2F2H(const double u, const double x1, const double x2,  const double yd) const;
	double SEXP(const double y, const double e) const;
	Eigen::MatrixXd source_matrix;
	Eigen::VectorXd rhs;
	Eigen::MatrixXd _src_matrix;
	void MakeMatrix(const double u, const double yd);
	void NewMakeMatrix(const double u);
	void MakeRhs(const double u);
	Eigen::MatrixXd MakeSrcMatrix() const;
	Eigen::MatrixXd if1_matrix, if2e_matrix, i1f2h_matrix, i2f2h_matrix;
	Eigen::VectorXd i1f2h_buf, if2e_buf;
	void fill_if1(const double u);
	void vect_if1_yd(const double u, const double yd, Eigen::VectorXd& buf);

	void fill_if2e(const double u); // OK
	void vect_if2e_yd(const double u, const double xd, const double yd, Eigen::VectorXd& buf);

	void fill_i1f2h(const double u);
	void vect_i1f2h(const double u, const double xd, Eigen::VectorXd& buf);
	void vect_i1f2h_yd(const double u, const double xd, const double yd, Eigen::VectorXd& buf);

	void fill_i2f2h(const double u);
	void vect_i2f2h_yd(const double u, const double yd, Eigen::VectorXd& buf); //TODO

	Eigen::VectorXd MakeGreenVector(const double u, const double xd, const double yd);

};
