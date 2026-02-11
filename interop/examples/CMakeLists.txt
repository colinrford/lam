
add_executable(finite_field_roots finite_field_roots.cpp)

target_link_libraries(finite_field_roots
  PRIVATE
    lam::interop
    lam::polynomial_nttp
    lam::ctbignum
)

target_compile_features(finite_field_roots PRIVATE cxx_std_23)
set_target_properties(finite_field_roots PROPERTIES CXX_SCAN_FOR_MODULES ON)

# Stress Test
add_executable(finite_field_stress_test finite_field_stress_test.cpp)

target_link_libraries(finite_field_stress_test
  PRIVATE
    lam::interop
    lam::polynomial_nttp
    lam::ctbignum
)

target_compile_features(finite_field_stress_test PRIVATE cxx_std_23)
set_target_properties(finite_field_stress_test PROPERTIES CXX_SCAN_FOR_MODULES ON)
