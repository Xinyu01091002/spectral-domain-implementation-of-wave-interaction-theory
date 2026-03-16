#pragma once

#include "mf12_cpp/spectral.hpp"

#include <filesystem>

namespace mf12_cpp {

void save_result_bundle(const std::filesystem::path& output_dir, const LoadedCase& loaded, const ResultBundle& result);

}  // namespace mf12_cpp
