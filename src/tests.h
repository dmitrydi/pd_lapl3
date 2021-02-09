/*
 * tests.h
 *
 *  Created on: 29 ���. 2021 �.
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

void TestNew_i12fh();
