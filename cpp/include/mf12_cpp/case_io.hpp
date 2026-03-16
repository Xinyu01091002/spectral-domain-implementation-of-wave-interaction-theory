#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace mf12_cpp {

struct Matrix {
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::vector<double> values;
};

struct CaseInputs {
  int order = 0;
  double g = 0.0;
  double h = 0.0;
  double Ux = 0.0;
  double Uy = 0.0;
  double Lx = 0.0;
  double Ly = 0.0;
  int Nx = 0;
  int Ny = 0;
  double t = 0.0;
  std::string subharmonic_mode;
};

struct CaseManifest {
  std::string case_id;
  std::string description;
  std::string purpose;
  CaseInputs inputs;
  std::map<std::string, std::filesystem::path> arrays;
  std::map<std::string, std::filesystem::path> reference_arrays;
  std::filesystem::path reference_metadata;
  std::map<std::string, double> tolerances;
};

struct LoadedCase {
  std::filesystem::path case_dir;
  CaseManifest manifest;
  std::map<std::string, Matrix> arrays;
  std::map<std::string, Matrix> reference_arrays;
};

CaseManifest load_manifest(const std::filesystem::path& case_dir);
LoadedCase load_case(const std::filesystem::path& case_dir);
Matrix load_csv_matrix(const std::filesystem::path& path);
bool validate_case(const LoadedCase& loaded, std::string* error_message);

}  // namespace mf12_cpp
