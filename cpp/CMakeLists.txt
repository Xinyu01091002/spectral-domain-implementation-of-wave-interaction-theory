cmake_minimum_required(VERSION 3.18)

project(mf12_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MF12_ENABLE_FFTW "Enable FFTW backend when a compatible FFTW3 install is available" OFF)
option(MF12_ENABLE_OPENMP "Enable OpenMP for hotspot parallelization when available" OFF)
set(MF12_FFTW_ROOT "C:/tools/fftw3-mingw-dll" CACHE PATH "Optional local FFTW installation root")
if(MF12_ENABLE_FFTW)
  find_package(FFTW3 CONFIG QUIET PATHS "${MF12_FFTW_ROOT}" NO_DEFAULT_PATH)
  if(NOT FFTW3_FOUND)
    find_path(
      MF12_FFTW_INCLUDE_DIR
      NAMES fftw3.h
      PATHS
        "${MF12_FFTW_ROOT}/include"
      NO_DEFAULT_PATH
    )
    find_library(
      MF12_FFTW_LIBRARY
      NAMES fftw3 libfftw3 libfftw3-3
      PATHS
        "${MF12_FFTW_ROOT}/lib"
      NO_DEFAULT_PATH
    )
    find_file(
      MF12_FFTW_DLL
      NAMES libfftw3-3.dll fftw3-3.dll
      PATHS
        "${MF12_FFTW_ROOT}/bin"
      NO_DEFAULT_PATH
    )
    if(MF12_FFTW_INCLUDE_DIR AND MF12_FFTW_LIBRARY)
      add_library(MF12::fftw3 UNKNOWN IMPORTED)
      set_target_properties(
        MF12::fftw3
        PROPERTIES
          IMPORTED_LOCATION "${MF12_FFTW_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${MF12_FFTW_INCLUDE_DIR}"
      )
      if(MF12_FFTW_DLL)
        set_target_properties(MF12::fftw3 PROPERTIES IMPORTED_RUNTIME_ARTIFACTS "${MF12_FFTW_DLL}")
      endif()
      set(FFTW3_FOUND TRUE)
    endif()
  endif()
endif()
if(MF12_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
endif()

add_executable(
  mf12_cpp
  src/case_io.cpp
  src/compare.cpp
  src/main.cpp
  src/result_io.cpp
  src/spectral.cpp
)

target_include_directories(mf12_cpp PRIVATE include)

if(MF12_ENABLE_FFTW AND FFTW3_FOUND)
  if(TARGET FFTW3::fftw3)
    target_link_libraries(mf12_cpp PRIVATE FFTW3::fftw3)
  else()
    target_link_libraries(mf12_cpp PRIVATE MF12::fftw3)
    target_include_directories(mf12_cpp PRIVATE "${MF12_FFTW_INCLUDE_DIR}")
  endif()
  target_compile_definitions(mf12_cpp PRIVATE MF12_HAVE_FFTW=1)
  if(MF12_FFTW_DLL)
    add_custom_command(
      TARGET mf12_cpp POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${MF12_FFTW_DLL}"
      "$<TARGET_FILE_DIR:mf12_cpp>"
    )
  endif()
endif()

if(MF12_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(mf12_cpp PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(mf12_cpp PRIVATE MF12_HAVE_OPENMP=1)
endif()
