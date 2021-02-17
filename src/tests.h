/*
 * tests.h
 *
 *  Created on: 29 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include "real_well.h"
#include "auxillary.h"
#include <map>
#include <utility>

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& s) {
  os << "{";
  bool first = true;
  for (const auto& x : s) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << x;
  }
  return os << "}";
}

template <class K, class V>
std::ostream& operator << (std::ostream& os, const std::map<K, V>& m) {
  os << "{" << std::endl;
  bool first = true;
  for (const auto& kv : m) {
    if (!first) {
      os << ", ";
      os << std::endl;
    }
    first = false;
    os << kv.first << ": " << kv.second;
  }
  return os <<std::endl << "}";
}

void SimpleTest();

void TestGringarten();

void TestCinco();

void TestCincoNew();

void TestNewSpeed();

void TestNew_i1f2h();

void TestNew_if2e();

void TestNew_if1();

void TestNew_i2f2h();

void TestGreen();

void GridPd();

void GridPdLapl();
