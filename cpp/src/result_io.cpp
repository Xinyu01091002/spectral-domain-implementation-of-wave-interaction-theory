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
  save_csv(output_dir / "u.csv", result.kinematics.u);
  save_csv(output_dir / "v.csv", result.kinematics.v);
  save_csv(output_dir / "w.csv", result.kinematics.w);
  save_csv(output_dir / "p.csv", result.kinematics.p);
  save_csv(output_dir / "phi_vol.csv", result.kinematics.phi_vol);
  save_csv(output_dir / "uV.csv", result.kinematics.uV);
  save_csv(output_dir / "vV.csv", result.kinematics.vV);
  save_csv(output_dir / "a_x.csv", result.kinematics.a_x);
  save_csv(output_dir / "a_y.csv", result.kinematics.a_y);

  std::ofstream out(output_dir / "result.json");
  if (!out) {
    throw std::runtime_error("Failed to write result.json");
  }
  out << std::setprecision(17);
#if defined(MF12_HAVE_FFTW)
  const char* fft_backend = "FFTW3 backward 2D complex DFT";
#else
  const char* fft_backend = "in-house radix-2 inverse FFT (fallback direct inverse DFT)";
#endif
  out << "{\n"
      << "  \"case_id\": \"" << loaded.manifest.case_id << "\",\n"
      << "  \"runtime\": {\n"
      << "    \"repeats\": " << result.runtime.repeats << ",\n"
      << "    \"mean_coefficient_s\": " << result.runtime.mean_coefficient_s << ",\n"
      << "    \"mean_linear_coefficient_s\": " << result.runtime.mean_linear_coefficient_s << ",\n"
      << "    \"mean_second_order_coefficient_s\": " << result.runtime.mean_second_order_coefficient_s << ",\n"
      << "    \"mean_third_order_coefficient_s\": " << result.runtime.mean_third_order_coefficient_s << ",\n"
      << "    \"mean_third_order_np2m_s\": " << result.runtime.mean_third_order_np2m_s << ",\n"
      << "    \"mean_third_order_2npm_s\": " << result.runtime.mean_third_order_2npm_s << ",\n"
      << "    \"mean_third_order_npmpp_s\": " << result.runtime.mean_third_order_npmpp_s << ",\n"
      << "    \"mean_reconstruction_s\": " << result.runtime.mean_reconstruction_s << ",\n"
      << "    \"mean_kinematics_s\": " << result.runtime.mean_kinematics_s << ",\n"
      << "    \"mean_total_s\": " << result.runtime.mean_total_s << ",\n"
      << "    \"best_total_s\": " << result.runtime.best_total_s << "\n"
      << "  },\n"
      << "  \"metadata\": {\n"
      << "    \"language\": \"cpp\",\n"
      << "    \"implementation\": \"mf12_cpp spectral-only order<=3 (third-order WIP)\",\n"
      << "    \"fft_backend\": \"" << fft_backend << "\"\n"
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
