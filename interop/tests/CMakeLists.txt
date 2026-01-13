
enable_testing()

add_executable(test_interop_roots roots_finite_field.cpp)

target_link_libraries(test_interop_roots 
  PRIVATE
    lam::interop
    lam::concepts
    lam::polynomial_nttp
    lam::ctbignum
)

target_compile_features(test_interop_roots PRIVATE cxx_std_23)
set_target_properties(test_interop_roots PROPERTIES CXX_SCAN_FOR_MODULES ON)

add_test(NAME Interop_Roots_Finite_Field COMMAND test_interop_roots)
