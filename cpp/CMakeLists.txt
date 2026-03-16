cmake_minimum_required(VERSION 3.18)

project(mf12_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(
  mf12_cpp
  src/case_io.cpp
  src/compare.cpp
  src/main.cpp
  src/result_io.cpp
  src/spectral.cpp
)

target_include_directories(mf12_cpp PRIVATE include)
