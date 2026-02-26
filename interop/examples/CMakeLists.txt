add_interop_executable(finite_field_roots)
add_interop_executable(finite_field_stress_test)

# ============================================================
# Compatibility Examples
# ============================================================

add_interop_executable(ctbignum_compatibility SOURCE compatibility/ctbignum_compat.cpp)
target_compile_options(ctbignum_compatibility PRIVATE -fconstexpr-steps=150000000 -fbracket-depth=1024)

add_interop_executable(roots_brute_force SOURCE compatibility/roots_brute_force.cpp)
add_interop_executable(root_finding_demo)
