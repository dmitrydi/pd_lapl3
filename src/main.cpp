//============================================================================
// Name        : pd_lapl3.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include "lapl_well.h"
#include "profile.h"
#include "real_well.h"
#include "auxillary.h"
#include <cmath>
#include <utility>
#include "tests.h"
using namespace std;

int main() {
	//SimpleTest(); OK
	// TestGringarten();
	// TestCinco(); OK
	//LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_, const double alpha_ = 0.);
	double xed = 10.;
	double xwd = 3.;
	double yed = 3;
	double ywd = yed/2.;
	double Fcd = 3.14;
	double u = 1.2;
	LaplWell lwell(xwd, xed, ywd, yed, Fcd);
	cout << lwell.make_if2e(u);
}
