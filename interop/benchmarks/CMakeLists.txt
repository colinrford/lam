# Interop Benchmarks
add_interop_executable(benchmark_ntt_pure)
add_interop_executable(benchmark_ntt_vs_naive)
add_interop_executable(benchmark_parallel_ntt)
add_interop_executable(benchmark_ntt_vector_vs_array)

# NTL Benchmark (if available)
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h PATHS /usr/local/include /opt/homebrew/include)
find_library(NTL_LIBRARY NAMES ntl PATHS /usr/local/lib /opt/homebrew/lib)

if(NTL_INCLUDE_DIR AND NTL_LIBRARY)
  message(STATUS "NTL found: enabling NTL benchmark in interop")
  add_interop_executable(benchmark_ntt_vs_ntl SOURCE benchmark_ntt_vs_ntl.cpp ntl_benchmark_impl.cpp LIBS ${NTL_LIBRARY})
  target_include_directories(benchmark_ntt_vs_ntl PRIVATE ${NTL_INCLUDE_DIR})
else()
  message(STATUS "NTL not found")
endif()
