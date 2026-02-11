/*
 *  finite_field_roots.cpp
 *    Demonstrates finding roots of polynomials over finite fields
 *    using lam::polynomial_nttp and lam::ctbignum
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop; // pulls in essential traits
import lam.interop.extra;

using namespace lam::cbn::literals;
using namespace lam::polynomial;

// Define a field GF(p) with p = 17 for this example
constexpr auto P = 17_Z;
using GF17 = decltype(lam::cbn::ZqElement<std::uint64_t, 17>(0));

void demo_quadratic()
{
    std::println("\n--- Solving x^2 - 2 = 0 in GF(17) ---");
    // Squares in GF(17):
    // 1^2=1, 2^2=4, 3^2=9, 4^2=16, 5^2=25=8, 6^2=36=2, 7^2=49=15, 8^2=64=13
    // So 6^2 = 2. Also (-6)^2 = 11^2 = 121 = 2.
    // Roots should be 6 and 11.

    // Polynomial: 1*x^2 + 0*x + (-2)
    // -2 mod 17 = 15
    constexpr GF17 c2 = 1_Z;
    constexpr GF17 c1 = 0_Z;
    constexpr GF17 c0 = 15_Z;

    // Using 'roots' function which dispatches to Berlekamp for finite fields
    univariate::polynomial_nttp<GF17, 2> poly{{c0, c1, c2}};
    
    lam::interop::print_poly(poly, "P(x)");

    auto results = lam::roots(poly);

    std::println("Found {} roots:", results.size());
    for (const auto& r : results)
    {
        std::println("  Root: {}, Multiplicity: {}", (std::uint64_t)r.value.data[0], r.multiplicity);
    }
}

void demo_cubic()
{
    std::println("\n--- Solving Cubic in GF(17) ---");
    // Let's construct a polynomial with known roots: (x-1)(x-2)(x-3)
    // = (x^2 - 3x + 2)(x - 3)
    // = x^3 - 3x^2 + 2x - 3x^2 + 9x - 6
    // = x^3 - 6x^2 + 11x - 6
    
    // In GF(17):
    // -6 = 11
    // 11 = 11
    // -6 = 11
    // So: x^3 + 11x^2 + 11x + 11
    
    GF17 one = 1_Z;
    GF17 eleven = 11_Z;
    
    univariate::polynomial_nttp<GF17, 3> poly{{eleven, eleven, eleven, one}};
     
    lam::interop::print_poly(poly, "Q(x)");
    
    auto results = lam::roots(poly);
    
    std::println("Found {} roots (expecting 1, 2, 3):", results.size());
    for (const auto& r : results)
    {
         std::println("  Root: {}, Multiplicity: {}", (std::uint64_t)r.value.data[0], r.multiplicity);
    }
}

void demo_quartic()
{
    std::println("\n--- Solving Quartic in GF(17) ---");
    // (x-1)(x-2)(x-3)(x-4)
    // Roots: 1, 2, 3, 4
    
    // Expand to get coeffs:
    // (x^2 - 3x + 2)(x^2 - 7x + 12)
    // x^4 - 7x^3 + 12x^2
    //     - 3x^3 + 21x^2 - 36x
    //            + 2x^2 - 14x + 24
    // ----------------------------
    // x^4 - 10x^3 + 35x^2 - 50x + 24
    
    // In GF(17):
    // -10 = 7
    // 35 = 1
    // -50 = -16 = 1
    // 24 = 7
    
    // So: x^4 + 7x^3 + x^2 + x + 7
    
    using GF = GF17;
    univariate::polynomial_nttp<GF, 4> poly{{
        GF(7_Z), // x^0
        GF(1_Z), // x^1
        GF(1_Z), // x^2
        GF(7_Z), // x^3
        GF(1_Z)  // x^4
    }};
    
    lam::interop::print_poly(poly, "R(x)");
    auto results = lam::roots(poly);
    
    std::println("Found {} roots (expecting 1, 2, 3, 4):", results.size());
    for (const auto& r : results)
    {
         std::println("  Root: {}, Multiplicity: {}", (std::uint64_t)r.value.data[0], r.multiplicity);
    }
}

void demo_mixed_degree()
{
    std::println("\n--- Solving Mixed Degree in GF(17) ---");
    // (x - 5)(x^2 - 3)
    // x^2 - 3 is irreducible in GF(17) because 3 is not a quadratic residue.
    // So we expect only ONE root: 5.
    
    // Expansion:
    // x(x^2 - 3) - 5(x^2 - 3)
    // x^3 - 3x - 5x^2 + 15
    // x^3 - 5x^2 - 3x + 15
    
    // In GF(17):
    // -5 = 12
    // -3 = 14
    // 15 = 15
    
    // So: x^3 + 12x^2 + 14x + 15
    
    using GF = GF17;
    univariate::polynomial_nttp<GF, 3> poly{{
        GF(15_Z), // x^0
        GF(14_Z), // x^1
        GF(12_Z), // x^2
        GF(1_Z)   // x^3
    }};
    
    lam::interop::print_poly(poly, "S(x)");
    auto results = lam::roots(poly);
    
    std::println("Found {} roots (expecting only 5):", results.size());
    for (const auto& r : results) 
    {
         std::println("  Root: {}, Multiplicity: {}", (std::uint64_t)r.value.data[0], r.multiplicity);
    }
}

int main()
{
    std::println("Finite Field Root Finding Demo");
    demo_quadratic();
    demo_cubic();
    demo_quartic();
    demo_mixed_degree();
    return 0;
}

