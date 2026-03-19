#pragma once

#include "mf12_cpp/spectral.hpp"

#include <map>
#include <string>

namespace mf12_cpp {

std::map<std::string, double> compare_result_to_reference(const Matrix& eta, const Matrix& phi, const LoadedCase& loaded);
std::map<std::string, double> compare_result_to_reference(const ResultBundle& result, const LoadedCase& loaded);
bool tolerances_pass(const std::map<std::string, double>& metrics, const std::map<std::string, double>& tolerances);

}  // namespace mf12_cpp
