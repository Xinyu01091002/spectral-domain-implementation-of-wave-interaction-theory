#pragma once

#include "mf12_cpp/case_io.hpp"

#include <map>
#include <string>

namespace mf12_cpp {

std::map<std::string, double> compare_result_to_reference(const Matrix& eta, const Matrix& phi, const LoadedCase& loaded);
bool tolerances_pass(const std::map<std::string, double>& metrics, const std::map<std::string, double>& tolerances);

}  // namespace mf12_cpp
