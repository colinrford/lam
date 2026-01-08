cmake_minimum_required(VERSION 3.31.6)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31.6)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

project(lam VERSION 0.1.260108 LANGUAGES CXX)

# Build individual libraries in dependency order
enable_testing()
set(LAM_CTBIGNUM_BuildTests ON CACHE BOOL "" FORCE)

add_subdirectory(concepts)
add_subdirectory(ctbignum)
add_subdirectory(polynomial_nttp)
add_subdirectory(bitvector)

# Build interop module (depends on all above)
add_subdirectory(interop)
