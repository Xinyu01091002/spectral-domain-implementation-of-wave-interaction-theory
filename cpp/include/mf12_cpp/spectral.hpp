#pragma once

#include "mf12_cpp/case_io.hpp"

#include <complex>
#include <map>
#include <string>
#include <vector>

namespace mf12_cpp {

struct RuntimeStats {
  int repeats = 1;
  bool warmup = false;
  double mean_coefficient_s = 0.0;
  double mean_linear_coefficient_s = 0.0;
  double mean_second_order_coefficient_s = 0.0;
  double mean_third_order_coefficient_s = 0.0;
  double mean_third_order_np2m_s = 0.0;
  double mean_third_order_2npm_s = 0.0;
  double mean_third_order_npmpp_s = 0.0;
  double mean_reconstruction_s = 0.0;
  double mean_kinematics_s = 0.0;
  double mean_total_s = 0.0;
  double best_total_s = 0.0;
};

struct KinematicsBundle {
  Matrix u;
  Matrix v;
  Matrix w;
  Matrix p;
  Matrix phi_vol;
  Matrix uV;
  Matrix vV;
  Matrix a_x;
  Matrix a_y;
};

struct ResultBundle {
  Matrix eta;
  Matrix phi;
  Matrix x;
  Matrix y;
  KinematicsBundle kinematics;
  RuntimeStats runtime;
  std::map<std::string, double> comparison;
};

ResultBundle run_case(const LoadedCase& loaded, int repeats = 1, bool warmup = false);
void dump_coefficients_for_case(const LoadedCase& loaded, const std::filesystem::path& output_dir);

}  // namespace mf12_cpp
