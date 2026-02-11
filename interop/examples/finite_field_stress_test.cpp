/*
 *  finite_field_stress_test.cpp
 *    Stress testing finite field root finding with corner cases
 *    and "interesting" polynomials.
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;
import lam.interop.extra;

using namespace lam::cbn::literals;
using namespace lam::polynomial;

// GF(17)
constexpr auto P = 17_Z;
using GF17 = decltype(lam::cbn::ZqElement<std::uint64_t, 17>(0));

void test_field_polynomial()
{
    std::println("\n--- Test: Field Polynomial x^17 - x in GF(17) ---");
    std::println("Expectation: All 17 elements of GF(17) are roots.");

    // x^17 - x
    univariate::polynomial_nttp<GF17, 17> poly;
    poly.coefficients[17] = GF17(1_Z);
    poly.coefficients[1] = GF17(16_Z); // -1 = 16
    
    // lam::interop::print_poly(poly, "P(x)"); // Too long to print nicely?

    auto results = lam::roots(poly);
    
    std::println("Found {} roots.", results.size());
    
    if (results.size() == 17) {
        std::println("[PASS] Found all 17 distinct roots.");
    } else {
        std::println("[FAIL] Expected 17 distinct roots, found {}.", results.size());
    }
}

void test_repeated_roots()
{
    std::println("\n--- Test: Repeated Roots (x-1)^17 = x^17 - 1 in GF(17) ---");
    std::println("Expectation: One root (1) with multiplicity 17.");
    std::println("Note: Derivative is identically zero here.");

    // x^17 - 1
    univariate::polynomial_nttp<GF17, 17> poly;
    poly.coefficients[17] = GF17(1_Z);
    poly.coefficients[0] = GF17(16_Z); // -1 = 16

    auto results = lam::roots(poly);
    
    std::println("Found {} unique roots.", results.size());
    for(const auto& r : results) {
        std::println("  Root: {}, Multiplicity: {}", (std::uint64_t)r.value.data[0], r.multiplicity);
    }

    if (results.size() == 1 && results[0].value == GF17(1_Z) && results[0].multiplicity == 17) {
        std::println("[PASS] Found root 1 with multiplicity 17.");
    } else if (results.size() > 0 && results[0].value == GF17(1_Z)) {
        std::println("[PARTIAL] Found root 1 but multiplicity {} (expected 17).", results[0].multiplicity);
    } else {
        std::println("[FAIL] Failed to find root 1.");
    }
}

void test_cyclotomic()
{
    std::println("\n--- Test: Cyclotomic Phi_8(x) = x^4 + 1 in GF(17) ---");
    std::println("Expectation: 4 distinct primitive 8th roots of unity (2, 8, 9, 15).");

    // x^4 + 1
    univariate::polynomial_nttp<GF17, 4> poly;
    poly.coefficients[4] = GF17(1_Z);
    poly.coefficients[0] = GF17(1_Z);
    
    lam::interop::print_poly(poly, "Phi_8(x)");

    auto results = lam::roots(poly);
    
    std::println("Found {} roots.", results.size());
    std::vector<std::uint64_t> expected = {2, 8, 9, 15};
    std::sort(expected.begin(), expected.end());
    
    std::vector<std::uint64_t> found;
    for(const auto& r : results) found.push_back((std::uint64_t)r.value.data[0]);
    std::sort(found.begin(), found.end());
    
    if (found == expected) {
        std::println("[PASS] Found exact set {{2, 8, 9, 15}}.");
    } else {
        std::print("[FAIL] Found: ");
        for(auto v : found) std::print("{} ", v);
        std::println("");
    }
}

int main()
{
    std::println("Finite Field Stress Test");
    test_field_polynomial();
    test_repeated_roots();
    test_cyclotomic();
    return 0;
}
