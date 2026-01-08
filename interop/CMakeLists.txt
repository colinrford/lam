add_library(lam_interop)

target_sources(lam_interop
  PUBLIC
    FILE_SET CXX_MODULES FILES
      src/lam-interop.cppm
)

target_compile_features(lam_interop PUBLIC cxx_std_23)

target_link_libraries(lam_interop
  PUBLIC
    lam::concepts
    lam::ctbignum
    lam::polynomial_nttp
    lam::bitvector
)

add_library(lam::interop ALIAS lam_interop)
