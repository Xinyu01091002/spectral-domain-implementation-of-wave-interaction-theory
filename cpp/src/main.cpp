#include "mf12_cpp/case_io.hpp"
#include "mf12_cpp/compare.hpp"
#include "mf12_cpp/result_io.hpp"
#include "mf12_cpp/spectral.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void print_usage() {
  std::cout
      << "mf12_cpp commands:\n"
      << "  inspect <case_dir>   Load a shared case and print a summary\n"
      << "  validate <case_dir>  Validate that a shared case matches the expected layout\n"
      << "  verify <case_dir> [output_dir] [order]  Run the C++ spectral path and compare to MATLAB reference\n"
      << "  dump-coeffs <case_dir> <output_dir> [order]  Dump coefficient arrays for debugging\n";
}

int run_inspect(const std::filesystem::path& case_dir) {
  const auto loaded = mf12_cpp::load_case(case_dir);
  const std::size_t n_components = loaded.arrays.at("a").values.size();
  const auto& inp = loaded.manifest.inputs;
  std::cout << "{\n"
            << "  \"case_id\": \"" << loaded.manifest.case_id << "\",\n"
            << "  \"description\": \"" << loaded.manifest.description << "\",\n"
            << "  \"purpose\": \"" << loaded.manifest.purpose << "\",\n"
            << "  \"grid\": {\"Nx\": " << inp.Nx << ", \"Ny\": " << inp.Ny << "},\n"
            << "  \"domain\": {\"Lx\": " << inp.Lx << ", \"Ly\": " << inp.Ly << "},\n"
            << "  \"order\": " << inp.order << ",\n"
            << "  \"component_count\": " << n_components << ",\n"
            << "  \"subharmonic_mode\": \"" << inp.subharmonic_mode << "\",\n"
            << "  \"has_reference\": true\n"
            << "}\n";
  return 0;
}

int run_validate(const std::filesystem::path& case_dir) {
  const auto loaded = mf12_cpp::load_case(case_dir);
  std::string message;
  const bool ok = mf12_cpp::validate_case(loaded, &message);
  if (!ok) {
    std::cerr << "Validation failed: " << message << "\n";
    return 1;
  }
  std::cout << "Validation passed for case '" << loaded.manifest.case_id << "'.\n";
  return 0;
}

int run_verify(const std::filesystem::path& case_dir, const std::filesystem::path& output_dir, int override_order) {
  auto loaded = mf12_cpp::load_case(case_dir);
  if (override_order > 0) {
    loaded.manifest.inputs.order = override_order;
  }
  auto result = mf12_cpp::run_case(loaded, 1, false);
  result.comparison = mf12_cpp::compare_result_to_reference(result.eta, result.phi, loaded);
  mf12_cpp::save_result_bundle(output_dir, loaded, result);
  const bool pass = mf12_cpp::tolerances_pass(result.comparison, loaded.manifest.tolerances);
  std::cout << "{\n"
            << "  \"case_id\": \"" << loaded.manifest.case_id << "\",\n"
            << "  \"output_dir\": \"" << output_dir.string() << "\",\n"
            << "  \"order\": " << loaded.manifest.inputs.order << ",\n"
            << "  \"pass\": " << (pass ? "true" : "false") << ",\n"
            << "  \"metrics\": {\n";
  bool first = true;
  for (const auto& [key, value] : result.comparison) {
    if (!first) {
      std::cout << ",\n";
    }
    std::cout << "    \"" << key << "\": " << std::setprecision(17) << value;
    first = false;
  }
  std::cout << "\n  }\n}\n";
  return pass ? 0 : 1;
}

int run_dump_coeffs(const std::filesystem::path& case_dir, const std::filesystem::path& output_dir, int override_order) {
  auto loaded = mf12_cpp::load_case(case_dir);
  if (override_order > 0) {
    loaded.manifest.inputs.order = override_order;
  }
  mf12_cpp::dump_coefficients_for_case(loaded, output_dir);
  std::cout << "{\n"
            << "  \"case_id\": \"" << loaded.manifest.case_id << "\",\n"
            << "  \"order\": " << loaded.manifest.inputs.order << ",\n"
            << "  \"output_dir\": \"" << output_dir.string() << "\"\n"
            << "}\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 3) {
      print_usage();
      return 1;
    }

    const std::string command = argv[1];
    const std::filesystem::path case_dir = argv[2];

    if (command == "inspect") {
      return run_inspect(case_dir);
    }
    if (command == "validate") {
      return run_validate(case_dir);
    }
    if (command == "verify") {
      std::filesystem::path output_dir;
      if (argc >= 4) {
        output_dir = argv[3];
      } else {
        output_dir = std::filesystem::path("outputs") / "cross_language_comparison" / "verify_cpp" / "minimal_run";
      }
      int override_order = 0;
      if (argc >= 5) {
        override_order = std::stoi(argv[4]);
      }
      return run_verify(case_dir, output_dir, override_order);
    }
    if (command == "dump-coeffs") {
      if (argc < 4) {
        print_usage();
        return 1;
      }
      const std::filesystem::path output_dir = argv[3];
      int override_order = 0;
      if (argc >= 5) {
        override_order = std::stoi(argv[4]);
      }
      return run_dump_coeffs(case_dir, output_dir, override_order);
    }

    print_usage();
    return 1;
  } catch (const std::exception& exc) {
    std::cerr << "mf12_cpp error: " << exc.what() << "\n";
    return 1;
  }
}
