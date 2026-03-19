#include "mf12_cpp/case_io.hpp"

#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace mf12_cpp {
namespace {

std::string read_text(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open file: " + path.string());
  }
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

std::string extract_section(const std::string& text, const std::string& key) {
  const std::string needle = "\"" + key + "\"";
  const auto key_pos = text.find(needle);
  if (key_pos == std::string::npos) {
    return {};
  }
  const auto brace_start = text.find('{', key_pos);
  if (brace_start == std::string::npos) {
    return {};
  }
  int depth = 0;
  for (std::size_t i = brace_start; i < text.size(); ++i) {
    if (text[i] == '{') {
      ++depth;
    } else if (text[i] == '}') {
      --depth;
      if (depth == 0) {
        return text.substr(brace_start, i - brace_start + 1);
      }
    }
  }
  return {};
}

std::string extract_string(const std::string& text, const std::string& key) {
  const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
  std::smatch match;
  if (std::regex_search(text, match, pattern)) {
    return match[1].str();
  }
  return {};
}

double extract_number(const std::string& text, const std::string& key) {
  const std::regex pattern("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
  std::smatch match;
  if (!std::regex_search(text, match, pattern)) {
    throw std::runtime_error("Missing numeric key: " + key);
  }
  return std::stod(match[1].str());
}

int extract_int(const std::string& text, const std::string& key) {
  return static_cast<int>(extract_number(text, key));
}

std::map<std::string, std::filesystem::path> extract_path_map(const std::string& text) {
  std::map<std::string, std::filesystem::path> out;
  const std::regex pattern("\"([^\"]+)\"\\s*:\\s*\"([^\"]+)\"");
  auto begin = std::sregex_iterator(text.begin(), text.end(), pattern);
  auto end = std::sregex_iterator();
  for (auto it = begin; it != end; ++it) {
    out[it->str(1)] = std::filesystem::path(it->str(2));
  }
  return out;
}

std::map<std::string, double> extract_tolerances(const std::string& text) {
  std::map<std::string, double> out;
  const std::regex pattern("\"([^\"]+)\"\\s*:\\s*([-+0-9.eE]+)");
  auto begin = std::sregex_iterator(text.begin(), text.end(), pattern);
  auto end = std::sregex_iterator();
  for (auto it = begin; it != end; ++it) {
    out[it->str(1)] = std::stod(it->str(2));
  }
  return out;
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> cells;
  std::stringstream ss(line);
  std::string cell;
  while (std::getline(ss, cell, ',')) {
    cells.push_back(cell);
  }
  if (!line.empty() && line.back() == ',') {
    cells.emplace_back();
  }
  return cells;
}

}  // namespace

CaseManifest load_manifest(const std::filesystem::path& case_dir) {
  const auto text = read_text(case_dir / "case.json");

  CaseManifest manifest;
  manifest.case_id = extract_string(text, "case_id");
  manifest.description = extract_string(text, "description");
  manifest.purpose = extract_string(text, "purpose");

  const auto inputs = extract_section(text, "inputs");
  manifest.inputs.order = extract_int(inputs, "order");
  manifest.inputs.g = extract_number(inputs, "g");
  manifest.inputs.h = extract_number(inputs, "h");
  manifest.inputs.Ux = extract_number(inputs, "Ux");
  manifest.inputs.Uy = extract_number(inputs, "Uy");
  manifest.inputs.Lx = extract_number(inputs, "Lx");
  manifest.inputs.Ly = extract_number(inputs, "Ly");
  manifest.inputs.Nx = extract_int(inputs, "Nx");
  manifest.inputs.Ny = extract_int(inputs, "Ny");
  manifest.inputs.t = extract_number(inputs, "t");
  try {
    manifest.inputs.z_kinematics = extract_number(inputs, "z_kinematics");
  } catch (const std::exception&) {
    manifest.inputs.z_kinematics = 0.0;
  }
  manifest.inputs.subharmonic_mode = extract_string(inputs, "subharmonic_mode");

  manifest.arrays = extract_path_map(extract_section(text, "arrays"));

  const auto ref = extract_section(text, "reference");
  manifest.reference_arrays = extract_path_map(extract_section(ref, "arrays"));
  manifest.reference_metadata = extract_string(ref, "metadata");

  manifest.tolerances = extract_tolerances(extract_section(text, "tolerances"));
  return manifest;
}

Matrix load_csv_matrix(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open CSV: " + path.string());
  }

  Matrix mat;
  std::string line;
  std::size_t expected_cols = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }
    const auto cells = split_csv_line(line);
    if (expected_cols == 0) {
      expected_cols = cells.size();
      mat.cols = expected_cols;
    } else if (cells.size() != expected_cols) {
      throw std::runtime_error("Inconsistent CSV column count in: " + path.string());
    }
    for (const auto& cell : cells) {
      mat.values.push_back(std::stod(cell));
    }
    ++mat.rows;
  }
  return mat;
}

LoadedCase load_case(const std::filesystem::path& case_dir) {
  LoadedCase loaded;
  loaded.case_dir = std::filesystem::absolute(case_dir);
  loaded.manifest = load_manifest(loaded.case_dir);

  for (const auto& [name, rel_path] : loaded.manifest.arrays) {
    loaded.arrays[name] = load_csv_matrix(loaded.case_dir / rel_path);
  }
  for (const auto& [name, rel_path] : loaded.manifest.reference_arrays) {
    loaded.reference_arrays[name] = load_csv_matrix(loaded.case_dir / rel_path);
  }
  return loaded;
}

bool validate_case(const LoadedCase& loaded, std::string* error_message) {
  auto fail = [&](const std::string& message) {
    if (error_message != nullptr) {
      *error_message = message;
    }
    return false;
  };

  const auto find_array = [&](const std::map<std::string, Matrix>& arrays, const std::string& name) -> const Matrix* {
    const auto it = arrays.find(name);
    return it == arrays.end() ? nullptr : &it->second;
  };

  const Matrix* a = find_array(loaded.arrays, "a");
  const Matrix* b = find_array(loaded.arrays, "b");
  const Matrix* kx = find_array(loaded.arrays, "kx");
  const Matrix* ky = find_array(loaded.arrays, "ky");
  if (a == nullptr || b == nullptr || kx == nullptr || ky == nullptr) {
    return fail("Case is missing one of a, b, kx, ky.");
  }

  const std::size_t n = a->values.size();
  if (b->values.size() != n || kx->values.size() != n || ky->values.size() != n) {
    return fail("Input component arrays do not have matching lengths.");
  }

  const Matrix* eta = find_array(loaded.reference_arrays, "eta");
  const Matrix* phi = find_array(loaded.reference_arrays, "phi");
  const Matrix* x = find_array(loaded.reference_arrays, "x");
  const Matrix* y = find_array(loaded.reference_arrays, "y");
  if (eta == nullptr || phi == nullptr || x == nullptr || y == nullptr) {
    return fail("Reference arrays are incomplete.");
  }

  if (eta->rows != static_cast<std::size_t>(loaded.manifest.inputs.Ny) ||
      eta->cols != static_cast<std::size_t>(loaded.manifest.inputs.Nx)) {
    return fail("eta dimensions do not match Ny x Nx from case inputs.");
  }
  if (phi->rows != eta->rows || phi->cols != eta->cols) {
    return fail("phi dimensions do not match eta.");
  }
  if (x->values.size() != static_cast<std::size_t>(loaded.manifest.inputs.Nx)) {
    return fail("x axis length does not match Nx.");
  }
  if (y->values.size() != static_cast<std::size_t>(loaded.manifest.inputs.Ny)) {
    return fail("y axis length does not match Ny.");
  }

  const char* kinematic_names[] = {"u", "v", "w", "p", "phi_vol", "uV", "vV", "a_x", "a_y"};
  for (const char* name : kinematic_names) {
    const Matrix* mat = find_array(loaded.reference_arrays, name);
    if (mat != nullptr) {
      if (mat->rows != eta->rows || mat->cols != eta->cols) {
        return fail(std::string(name) + " dimensions do not match eta.");
      }
    }
  }
  return true;
}

}  // namespace mf12_cpp
