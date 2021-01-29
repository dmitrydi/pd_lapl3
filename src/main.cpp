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
using namespace std;

int main() {
	double xed = 10;
	double xwd = 0.5*xed;
	double yed = 10;
	double ywd = 0.5*yed;
	double Fcd = 10.;
	RealWell well(xwd, xed, ywd, yed, Fcd);
	vector<double> tds = LogSpaced(0.0001, 1000., 60);
	for (auto td: tds) {
		cout << td << "," << well.pwd(td) << endl;
		//break;
	}

}
