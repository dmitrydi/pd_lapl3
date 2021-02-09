/*
 * tests.cpp
 *
 *  Created on: 29 џэт. 2021 у.
 *      Author: Dmitry_Di
 */
#include "tests.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

using namespace std;

void SimpleTest() {
	double xed = 10;
	double xwd = 0.5*xed;
	double yed = 10;
	double ywd = 0.5*yed;
	double Fcd = 10.;
	RealWell well(xwd, xed, ywd, yed, Fcd);
	vector<double> tds = LogSpaced(0.0001, 1000., 60);
	LOG_DURATION("total");
	for (auto td: tds) {
		cout << td << "," << well.pwd(td) << endl;
	}
}

pair<vector<double>, map<int, vector<double>>> ParseGringarten(const string& datafile) {
	ifstream input(datafile);
	map<int, vector<double>> ref_data;
	vector<double> tds;
	if (input) {
		string line;
		getline(input, line); // xe/xf line
		stringstream is(line);
		string s;
		getline(is, s, ';'); // dummy getine;

		vector<int> xe_xf;

		while (getline(is, s, ';')) { // initializing xe/xf data
			xe_xf.push_back(stoi(s));
		}
		getline(input, line); // dummy getline
		while (getline(input, line)) {
			istringstream sline(line);
			getline(sline, s, ';');
			tds.push_back(stod(s));
			for (int xexf : xe_xf) {
				getline(sline, s, ';');
				ref_data[xexf].push_back(stod(s));
			}
		}
	}
	return {tds, ref_data};
}

void TestGringarten() {
	const string datafile = "./test_data/gringarten_rectangular_reservoir_inf_cond_c.csv";
	auto ref_data = ParseGringarten(datafile);
	vector<double> tdxes = ref_data.first;
	map<int, vector<double>> data = ref_data.second;
	map<int, vector<double>> calculated;
	const double Fcd = 1e+7;
	for (auto& [ixe, pd] : data) {
		double xed = 2.*static_cast<double>(ixe);
		double yed = xed;
		double xwd = xed*0.5;
		double ywd = 0.5*yed;
		RealWell well(xwd, xed, ywd, yed, Fcd);
		for (size_t i = 0; i < tdxes.size(); ++i) {
			double tdxe = tdxes[i];
			double tdxf = static_cast<double>(xed*xed)*tdxe;
			double ans = well.pwd(tdxf);
			calculated[ixe].push_back(ans);
			double eps = abs(ans-pd[i])/pd[i];
			cout << "xe/xf: " << ixe << " td_xf: " << tdxf << " calculated: " << ans << " expected: " << pd[i] << " eps: " << eps << endl;
			//break;
		}
		break;
	}
	ofstream fout("./test_data/test_gringarten_out.csv");
	for (size_t i = 0; i < tdxes.size(); ++i) {
		for (auto it = calculated.begin(); it != calculated.end(); ++it) {
			fout << calculated[it->first][i] << ';';
		}
		fout << '\n';
	}
}

void TestCinco() {
	const double PI = 3.141592653589793;
	const string datafile = "./test_data/cinco_infinite_reservoir.csv";
	const vector<double> Fcds = {PI*0.2, PI, 2.*PI, 10.*PI, 20.*PI, 100.*PI};
	ifstream dfile(datafile);
	vector<double> tds;
	vector<vector<double>> pwd_expected(Fcds.size());
	vector<vector<double>> calculated(Fcds.size());
	if (dfile) {
		string line;
		while (getline(dfile, line)) {
			if (line[0] == '/') continue;
			istringstream is(line);
			string s;
			getline(is, s, ';');
			tds.push_back(stod(s));
			for (size_t i = 0; i < Fcds.size(); ++i) {
				getline(is, s, ';');
				pwd_expected[i].push_back(stod(s));
			}
		}
	}
	const double xed = 100.;
	const double yed = 100.;
	const double xwd = xed*0.5;
	const double ywd = yed*0.5;
	for (size_t i = 0; i < Fcds.size(); ++i) {
		double Fcd = Fcds[i];
		RealWell well(xwd, xed, ywd, yed, Fcd);
		for (size_t j = 0; j < tds.size(); ++j) {
			double td  = tds[j];
			double pwd = well.pwd(td);
			double pwd_e = pwd_expected[i][j];
			double eps = abs(pwd - pwd_e)/pwd_e;
			cout << "Fcd: " << Fcd << " td: " << td << " pwd: " << pwd << " expected: " << pwd_e << " eps: " << eps << endl;
			calculated[i].push_back(pwd);
		}
	}

	ofstream out("./test_data/cinco_test_out.csv");
	for (size_t i = 0; i < tds.size(); ++i) {
		out << tds[i] << ';';
		for (size_t j = 0; j < Fcds.size(); ++j) {
			out << calculated[j][i];
			if (j < Fcds.size() - 1) out << ';';
		}
		out << '\n';
	}
}

void TestNew_i1f2h() {
	double xed = 10.;
	double xwd = 3.;
	double yed = 3;
	double ywd = yed/2.;
	double Fcd = 3.14;
	double u = 1.2;
	LaplWell lwell(xwd, xed, ywd, yed, Fcd);
	Eigen::MatrixXd old_ans(2*NSEG, 2*NSEG), new_ans(2*NSEG, 2*NSEG);
	int N = 100;
	{
		LOG_DURATION("old");
		for (int i = 0; i<N; ++i)
			old_ans = lwell.make_i1f2h( u);
	}
	{
		LOG_DURATION("new");
		for (int i = 0; i<N; ++i)
			new_ans = lwell.show_i1f2h_matrix(u);
	}
	cout << old_ans << endl;
	cout << "----------------\n";
	cout << new_ans << endl;
	Eigen::MatrixXd delta = old_ans-new_ans;
	cout << "----------------\n";
	double max_delta = 0;
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			if (abs(delta(i,j))>max_delta) max_delta = abs(delta(i,j));
		}
	}
	cout << max_delta;
}

