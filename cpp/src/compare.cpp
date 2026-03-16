#include "mf12_cpp/compare.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace mf12_cpp {
namespace {

double vector_norm(const std::vector<double>& values) {
  double sum = 0.0;
  for (double value : values) {
    sum += value * value;
  }
  return std::sqrt(sum);
}

std::map<std::string, double> compare_fields(const Matrix& candidate, const Matrix& reference, const std::string& prefix) {
  if (candidate.rows != reference.rows || candidate.cols != reference.cols) {
    throw std::runtime_error("Field dimensions do not match for comparison.");
  }
  std::vector<double> diff(candidate.values.size(), 0.0);
  double max_abs = 0.0;
  double sum_sq = 0.0;
  for (std::size_t i = 0; i < candidate.values.size(); ++i) {
    diff[i] = candidate.values[i] - reference.values[i];
    const double abs_val = std::abs(diff[i]);
    if (abs_val > max_abs) {
      max_abs = abs_val;
    }
    sum_sq += diff[i] * diff[i];
  }
  const double denom = std::max(vector_norm(reference.values), std::numeric_limits<double>::epsilon());
  return {
      {prefix + "_max_abs_err", max_abs},
      {prefix + "_rms_err", std::sqrt(sum_sq / static_cast<double>(diff.size()))},
      {prefix + "_relative_l2_err", vector_norm(diff) / denom},
  };
}

}  // namespace

std::map<std::string, double> compare_result_to_reference(const Matrix& eta, const Matrix& phi, const LoadedCase& loaded) {
  std::map<std::string, double> metrics;
  const Matrix& eta_ref = loaded.reference_arrays.at("eta");
  const Matrix& phi_ref = loaded.reference_arrays.at("phi");
  const auto eta_metrics = compare_fields(eta, eta_ref, "eta");
  const auto phi_metrics = compare_fields(phi, phi_ref, "phi");
  metrics.insert(eta_metrics.begin(), eta_metrics.end());
  metrics.insert(phi_metrics.begin(), phi_metrics.end());
  return metrics;
}

bool tolerances_pass(const std::map<std::string, double>& metrics, const std::map<std::string, double>& tolerances) {
  for (const auto& [key, limit] : tolerances) {
    const auto it = metrics.find(key);
    if (it != metrics.end() && it->second > limit) {
      return false;
    }
  }
  return true;
}

}  // namespace mf12_cpp
