#include "mf12_cpp/result_io.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace mf12_cpp {
namespace {

void save_csv(const std::filesystem::path& path, const Matrix& mat) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Failed to write CSV: " + path.string());
  }
  out << std::setprecision(17);
  for (std::size_t r = 0; r < mat.rows; ++r) {
    for (std::size_t c = 0; c < mat.cols; ++c) {
      if (c > 0) {
        out << ',';
      }
      out << mat.values[r * mat.cols + c];
    }
    out << '\n';
  }
}

}  // namespace

void save_result_bundle(const std::filesystem::path& output_dir, const LoadedCase& loaded, const ResultBundle& result) {
  std::filesystem::create_directories(output_dir);
  save_csv(output_dir / "eta.csv", result.eta);
  save_csv(output_dir / "phi.csv", result.phi);
  save_csv(output_dir / "x.csv", result.x);
  save_csv(output_dir / "y.csv", result.y);

  std::ofstream out(output_dir / "result.json");
  if (!out) {
    throw std::runtime_error("Failed to write result.json");
  }
  out << std::setprecision(17);
  out << "{\n"
      << "  \"case_id\": \"" << loaded.manifest.case_id << "\",\n"
      << "  \"runtime\": {\n"
      << "    \"repeats\": " << result.runtime.repeats << ",\n"
      << "    \"mean_coefficient_s\": " << result.runtime.mean_coefficient_s << ",\n"
      << "    \"mean_reconstruction_s\": " << result.runtime.mean_reconstruction_s << ",\n"
      << "    \"mean_total_s\": " << result.runtime.mean_total_s << ",\n"
      << "    \"best_total_s\": " << result.runtime.best_total_s << "\n"
      << "  },\n"
      << "  \"metadata\": {\n"
      << "    \"language\": \"cpp\",\n"
      << "    \"implementation\": \"mf12_cpp spectral-only order<=3 (third-order WIP)\",\n"
      << "    \"fft_backend\": \"built-in separable inverse DFT\"\n"
      << "  },\n"
      << "  \"comparison\": {\n";
  bool first = true;
  for (const auto& [key, value] : result.comparison) {
    if (!first) {
      out << ",\n";
    }
    out << "    \"" << key << "\": " << value;
    first = false;
  }
  out << "\n  }\n"
      << "}\n";
}

}  // namespace mf12_cpp
