/*
 * auxillary.h
 *
 *  Created on: 7 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <cmath>
#include <vector>

template <typename T>
std::vector<T> LogSpaced(const T min, const T max, const int steps) {
	std::vector<double> ans(steps);
	double ln_min = std::log(min);
	double ln_max = std::log(max);
	double d = (ln_max - ln_min)/(steps - 1);
	for (int i = 0; i < steps; ++i) {
		ans[i] = std::exp(ln_min + i*d);
	}
	return ans;
}

template <typename T>
std::vector<T> LinSpaced(const T min, const T max, const int steps) {
	std::vector<T> ans(steps);
	double d = (max - min)/(steps - 1);
	for (int i = 0; i < steps; ++i) {
		ans[i] = min + i*d;
	}
	return ans;
}
