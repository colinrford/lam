/*
 *  poly_rem_finite_field.cpp – Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit tests for poly_rem over finite field coefficient types (ZqElement)
 *
 *  Tests are structured as consteval functions verified by static_assert
 *  (compile-time) and also exercised with volatile-sourced dividend
 *  coefficients to force runtime execution (runtime coverage).
 *
 *  Tests resistant to compile-time evaluation: none identified; all ZqElement
 *    arithmetic (including /) is constexpr throughout.
 *  Tests resistant to runtime: none — Mod is always a NTTP (compile-time),
 *    but dividend coefficients can be volatile-sourced for runtime paths.
 *
 *  Fields and moduli exercised:
 *    GF(7)   — modulus x²+1           (odd prime, monic, degree 2)
 *    GF(5)   — modulus x²+2           (different odd prime, monic, degree 2)
 *    GF(7)   — modulus x³+x+1         (degree 3, no roots in GF(7))
 *    GF(2)   — modulus x²+x+1         (characteristic 2; unique algebraic behaviour)
 *    GF(7)   — modulus 3x²+1          (non-monic: exercises lead-coefficient division)
 *    GF(17)  — modulus x²+3           (larger prime; x²+3 irreducible over GF(17))
 *    GF(7)   — modulus x⁴+1           (degree-4, exercises more loop iterations)
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam;
using namespace lam::cbn;

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)

using GF2 = ZqElement<std::uint64_t, 2>;
using GF5 = ZqElement<std::uint64_t, 5>;
using GF7 = ZqElement<std::uint64_t, 7>;
using GF17 = ZqElement<std::uint64_t, 17>;

// ============================================================
// GF(7) — modulus x²+1
// (no root in GF(7): k²+1 ≠ 0 for all k in {0,...,6})
// ============================================================

consteval bool test_poly_rem_gf7_m_lt_n_constant()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 0> a{{three}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == three && r[1] == zero;
}

consteval bool test_poly_rem_gf7_m_lt_n_linear()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 1> a{{two, three}}; // 3x + 2
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == two && r[1] == three;
}

consteval bool test_poly_rem_gf7_exact_linear()
{
  // x²-1 = x²+6 in GF(7) = (x+1)(x-1); at x=-1: 0
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 1> mod{{one, one}};
  constexpr polynomial_nttp<GF7, 2> a{{six, zero, one}};
  constexpr auto r = poly_rem<GF7, 1, mod>(a);
  return r[0] == zero;
}

consteval bool test_poly_rem_gf7_quadratic_reduces()
{
  // x²+x+1 = 1·(x²+1) + x → remainder x → {0,1}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 2> a{{one, one, one}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == zero && r[1] == one;
}

consteval bool test_poly_rem_gf7_cubic_quadratic()
{
  // x³+2x+3 mod (x²+1): x·(x²+1)=x³+x; subtract: (2x+3)-x = x+3 → {3,1}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 3> a{{three, two, zero, one}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == three && r[1] == one;
}

consteval bool test_poly_rem_gf7_quartic_quadratic()
{
  // x⁴+x²+1 = x²·(x²+1) + 1 → remainder 1 → {1,0}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 4> a{{one, zero, one, zero, one}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_rem_gf7_quartic_linear_remainder()
{
  // x⁴+x mod (x²+1) = x+1 → {1,1}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 4> a{{zero, one, zero, zero, one}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == one && r[1] == one;
}

consteval bool test_poly_rem_gf7_zero_dividend()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 3> a{};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == zero && r[1] == zero;
}

// ============================================================
// GF(5) — modulus x²+2
// (no root in GF(5): k²+2 ≠ 0 for all k in {0,...,4})
// ============================================================

consteval bool test_poly_rem_gf5_constant_passthrough()
{
  constexpr GF5 zero(0);
  constexpr GF5 two(2);
  constexpr GF5 four(4);
  constexpr polynomial_nttp<GF5, 2> mod{{two, zero, GF5(1)}};
  constexpr polynomial_nttp<GF5, 0> a{{four}};
  constexpr auto r = poly_rem<GF5, 2, mod>(a);
  return r[0] == four && r[1] == zero;
}

consteval bool test_poly_rem_gf5_quadratic_reduces()
{
  // x²+x+1 = 1·(x²+2) + x-1 = 1·(x²+2) + x+4 → {4,1}
  constexpr GF5 zero(0);
  constexpr GF5 one(1);
  constexpr GF5 four(4);
  constexpr polynomial_nttp<GF5, 2> mod{{GF5(2), zero, one}};
  constexpr polynomial_nttp<GF5, 2> a{{one, one, one}};
  constexpr auto r = poly_rem<GF5, 2, mod>(a);
  return r[0] == four && r[1] == one;
}

consteval bool test_poly_rem_gf5_cubic_quadratic()
{
  // 2x³+x²+3x+4 mod (x²+2): quotient 2x+1, remainder 4x+2 → {2,4}
  constexpr GF5 zero(0);
  constexpr GF5 one(1);
  constexpr GF5 two(2);
  constexpr GF5 three(3);
  constexpr GF5 four(4);
  constexpr polynomial_nttp<GF5, 2> mod{{two, zero, one}};
  constexpr polynomial_nttp<GF5, 3> a{{four, three, one, two}};
  constexpr auto r = poly_rem<GF5, 2, mod>(a);
  return r[0] == two && r[1] == four;
}

// ============================================================
// GF(7) — cubic modulus x³+x+1
// (no root in GF(7): k³+k+1 ≠ 0 for all k)
// ============================================================

consteval bool test_poly_rem_gf7_cubic_mod_m_lt_n()
{
  // linear dividend mod deg-3 poly: passthrough to degree-2 result type
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 five(5);
  constexpr polynomial_nttp<GF7, 3> mod{{one, one, zero, one}};
  constexpr polynomial_nttp<GF7, 1> a{{five, two}};
  constexpr auto r = poly_rem<GF7, 3, mod>(a);
  return r[0] == five && r[1] == two && r[2] == zero;
}

consteval bool test_poly_rem_gf7_cubic_mod_degree3()
{
  // x³+1 mod (x³+x+1): 1·(x³+x+1) + (-x) → 6x → {0,6,0}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 3> mod{{one, one, zero, one}};
  constexpr polynomial_nttp<GF7, 3> a{{one, zero, zero, one}};
  constexpr auto r = poly_rem<GF7, 3, mod>(a);
  return r[0] == zero && r[1] == six && r[2] == zero;
}

// ============================================================
// GF(2) — characteristic 2, modulus x²+x+1
//
// Irreducibility: 0²+0+1=1≠0, 1²+1+1=3≡1≠0 in GF(2). ✓
//
// Characteristic 2 properties exercised here:
//   -x = x for all x (subtraction = addition)
//   1+1 = 0
//   Leading coefficient is always 1, so division by lead is trivial
//
// x²+x+1 is the minimal polynomial of the primitive element of GF(4).
//
// Expected remainders:
//   x³ mod (x²+x+1): (x+1)(x²+x+1)+1 = x³+1, so remainder = 1
//   x⁴ mod (x²+x+1): (x²+x)(x²+x+1)+x, so remainder = x
//   x²+x mod (x²+x+1): 1·(x²+x+1)+(−1) = 1·(x²+x+1)+1, remainder = 1
// ============================================================

consteval bool test_poly_rem_gf2_constant_passthrough()
{
  // constant 1 mod (x²+x+1): deg(a)=0 < N=2, passthrough
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 0> a{{one}};
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_rem_gf2_linear_passthrough()
{
  // x+1 mod (x²+x+1): deg(a)=1 < N=2, passthrough
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 1> a{{one, one}};
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == one && r[1] == one;
}

consteval bool test_poly_rem_gf2_x2_plus_x()
{
  // x²+x mod (x²+x+1) = 1 (since x²+x = (x²+x+1) + 1 in GF(2))
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 2> a{{zero, one, one}}; // x²+x
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_rem_gf2_cubic()
{
  // x³ mod (x²+x+1) = 1
  // (x+1)(x²+x+1) = x³+x²+x+x²+x+1 = x³+1 in GF(2); so x³ = (x+1)(x²+x+1)+1
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 3> a{{zero, zero, zero, one}}; // x³
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_rem_gf2_quartic()
{
  // x⁴ mod (x²+x+1) = x
  // x⁴ = (x²+x)(x²+x+1) + x; verify: (x²+x)(x²+x+1) = x⁴+x³+x²+x³+x²+x = x⁴+x in GF(2)
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 4> a{{zero, zero, zero, zero, one}}; // x⁴
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == zero && r[1] == one;
}

consteval bool test_poly_rem_gf2_reconstruction()
{
  // x³ = (x+1)·(x²+x+1) + 1; verify q*mod + r = a
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 2> mod{{one, one, one}};
  constexpr polynomial_nttp<GF2, 3> a{{zero, zero, zero, one}};
  constexpr auto r = poly_rem<GF2, 2, mod>(a);
  // quotient = x+1
  constexpr polynomial_nttp<GF2, 1> q{{one, one}};
  constexpr auto reconstructed = q * mod + r;
  return reconstructed[0] == a[0] && reconstructed[1] == a[1] && reconstructed[2] == a[2] && reconstructed[3] == a[3];
}

// ============================================================
// Non-monic modulus — GF(7), modulus 3x²+1
//
// The lead coefficient is 3 (not 1). poly_rem divides by lead via
//   R factor = rem[i] / lead
// using ZqElement's operator/, exercising a code path untested by
// all monic moduli. 3⁻¹ = 5 in GF(7) (3·5=15≡1).
//
// Expected remainders:
//   (x³+x+1) mod (3x²+1):
//     i=3: factor=1/3=5; rem[1] -= 5*1 → 1-5=-4=3; rem[3] -= 5*3=15=1 → 0
//     remainder = 3x+1 → {1,3}
//   (3x²+x+1) mod (3x²+1):
//     i=2: factor=3/3=1; rem[0] -= 1*1=0; rem[1]-=0; rem[2]-=3=0
//     remainder = x → {0,1}
// ============================================================

consteval bool test_poly_rem_gf7_nonmonic_cubic()
{
  // (x³+x+1) mod (3x²+1) in GF(7) = 3x+1 → {1,3}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, three}};  // 3x²+1
  constexpr polynomial_nttp<GF7, 3> a{{one, one, zero, one}}; // x³+x+1
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == one && r[1] == three;
}

consteval bool test_poly_rem_gf7_nonmonic_m_equals_n()
{
  // (3x²+x+1) mod (3x²+1) in GF(7): quotient=1, remainder=x → {0,1}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, three}}; // 3x²+1
  constexpr polynomial_nttp<GF7, 2> a{{one, one, three}};    // 3x²+x+1
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == zero && r[1] == one;
}

consteval bool test_poly_rem_gf7_nonmonic_exact()
{
  // (3x²+1)·(x+2) = 3x³+6x²+x+2 in GF(7) — exact divisibility
  // remainder should be zero → {0,0}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, three}};   // 3x²+1
  constexpr polynomial_nttp<GF7, 3> a{{two, one, six, three}}; // 3x³+6x²+x+2
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == zero && r[1] == zero;
}

// ============================================================
// GF(17) — modulus x²+3
//
// Irreducibility over GF(17): need -3 to be a quadratic non-residue.
// -3 ≡ 14 (mod 17); 14^8 ≡ (-3)^8 = 3^8.
// 3^2=9, 3^4=81≡13, 3^8=13^2=169≡16≡-1 (mod 17).
// Since 14^((17-1)/2) = -1 ≢ 1, 14 is a QNR mod 17 → x²+3 irreducible ✓
//
// Expected remainders:
//   (x³+x+1) mod (x²+3):
//     i=3: factor=1; rem[1] -= 3 → 1-3=-2=15; remainder = 15x+1
//   (x⁴+2x+5) mod (x²+3):
//     i=4: factor=1; rem[2] -= 3 → 0-3=14; i=2: factor=14; rem[0] -= 14*3=42=8 → 5-8=-3=14
//     remainder = 2x+14
// ============================================================

consteval bool test_poly_rem_gf17_constant_passthrough()
{
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 three(3);
  constexpr GF17 seven(7);
  constexpr polynomial_nttp<GF17, 2> mod{{three, zero, one}}; // x²+3
  constexpr polynomial_nttp<GF17, 0> a{{seven}};
  constexpr auto r = poly_rem<GF17, 2, mod>(a);
  return r[0] == seven && r[1] == zero;
}

consteval bool test_poly_rem_gf17_cubic()
{
  // (x³+x+1) mod (x²+3) in GF(17): remainder = 15x+1 → {1,15}
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 three(3);
  constexpr GF17 fifteen(15);
  constexpr polynomial_nttp<GF17, 2> mod{{three, zero, one}};
  constexpr polynomial_nttp<GF17, 3> a{{one, one, zero, one}};
  constexpr auto r = poly_rem<GF17, 2, mod>(a);
  return r[0] == one && r[1] == fifteen;
}

consteval bool test_poly_rem_gf17_quartic()
{
  // (x⁴+2x+5) mod (x²+3) in GF(17): remainder = 2x+14 → {14,2}
  // Step 1 (i=4): factor=1; rem[2] -= 3 → 14; i=3: factor=0; i=2: factor=14; rem[0] -= 42=8 → 5-8=14
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 two(2);
  constexpr GF17 three(3);
  constexpr GF17 five(5);
  constexpr GF17 fourteen(14);
  constexpr polynomial_nttp<GF17, 2> mod{{three, zero, one}};
  constexpr polynomial_nttp<GF17, 4> a{{five, two, zero, zero, one}};
  constexpr auto r = poly_rem<GF17, 2, mod>(a);
  return r[0] == fourteen && r[1] == two;
}

consteval bool test_poly_rem_gf17_reconstruction()
{
  // x³+x+1 = x·(x²+3) + (15x+1); verify
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 three(3);
  constexpr GF17 fifteen(15);
  constexpr polynomial_nttp<GF17, 2> mod{{three, zero, one}};
  constexpr polynomial_nttp<GF17, 3> a{{one, one, zero, one}};
  constexpr auto r = poly_rem<GF17, 2, mod>(a);      // {1,15}
  constexpr polynomial_nttp<GF17, 1> q{{zero, one}}; // x
  constexpr auto reconstructed = q * mod + r;
  return reconstructed[0] == a[0] && reconstructed[1] == a[1] && reconstructed[2] == a[2] && reconstructed[3] == a[3];
}

// ============================================================
// Degree-4 modulus — GF(7), modulus x⁴+1
//
// x⁴+1 may not be irreducible over GF(7), but poly_rem requires only
// that Mod is a nonzero polynomial — irreducibility is not needed.
//
// Expected remainders (using x⁴ ≡ -1 = 6 under this modulus):
//   x⁵+2x+3 mod (x⁴+1):
//     x⁵ = x·x⁴ ≡ -x; so x⁵+2x+3 ≡ -x+2x+3 = x+3 → {3,1,0,0}
//   x⁷+1 mod (x⁴+1):
//     x⁴≡-1, x⁵≡-x, x⁶≡-x², x⁷≡-x³=6x³; so x⁷+1 ≡ 6x³+1 → {1,0,0,6}
//   x⁴+3x+2 mod (x⁴+1): remainder = (3x+2)+(−1) = 3x+1 → {1,3,0,0}
// ============================================================

consteval bool test_poly_rem_gf7_degree4_mod_quintic()
{
  // x⁵+2x+3 mod (x⁴+1) = x+3 → {3,1,0,0}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 4> mod{{one, zero, zero, zero, one}}; // x⁴+1
  constexpr polynomial_nttp<GF7, 5> a{{three, two, zero, zero, zero, one}};
  constexpr auto r = poly_rem<GF7, 4, mod>(a);
  return r[0] == three && r[1] == one && r[2] == zero && r[3] == zero;
}

consteval bool test_poly_rem_gf7_degree4_mod_x7()
{
  // x⁷+1 mod (x⁴+1) in GF(7):
  // x⁷ = x³·x⁴ ≡ x³·(-1) = -x³ = 6x³; remainder = 6x³+1 → {1,0,0,6}
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 4> mod{{one, zero, zero, zero, one}};
  constexpr polynomial_nttp<GF7, 7> a{{one, zero, zero, zero, zero, zero, zero, one}};
  constexpr auto r = poly_rem<GF7, 4, mod>(a);
  return r[0] == one && r[1] == zero && r[2] == zero && r[3] == six;
}

consteval bool test_poly_rem_gf7_degree4_mod_exact_multiple()
{
  // (x⁴+1)·(x²+3) = x⁶+3x⁴+x²+3 in GF(7); remainder = 0
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 4> mod{{one, zero, zero, zero, one}};
  // (x⁴+1)(x²+3) = x⁶+3x⁴+x²+3 → {3,0,1,0,3,0,1}
  constexpr polynomial_nttp<GF7, 6> a{{three, zero, one, zero, three, zero, one}};
  constexpr auto r = poly_rem<GF7, 4, mod>(a);
  return r[0] == zero && r[1] == zero && r[2] == zero && r[3] == zero;
}

// ============================================================
// Reconstruction: GF(7) x²+1 and GF(5) x²+2
// ============================================================

consteval bool test_poly_rem_gf7_reconstruction()
{
  // x³+2x+3 = x·(x²+1) + (x+3)
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr polynomial_nttp<GF7, 2> mod{{one, zero, one}};
  constexpr polynomial_nttp<GF7, 3> a{{three, two, zero, one}};
  constexpr auto r = poly_rem<GF7, 2, mod>(a);
  constexpr polynomial_nttp<GF7, 1> q{{zero, one}};
  constexpr auto reconstructed = (q * mod) + r;
  return reconstructed[0] == a[0] && reconstructed[1] == a[1] && reconstructed[2] == a[2] && reconstructed[3] == a[3];
}

consteval bool test_poly_rem_gf5_reconstruction()
{
  // 2x³+x²+3x+4 = (2x+1)·(x²+2) + (4x+2)
  constexpr GF5 zero(0);
  constexpr GF5 one(1);
  constexpr GF5 two(2);
  constexpr GF5 three(3);
  constexpr GF5 four(4);
  constexpr polynomial_nttp<GF5, 2> mod{{two, zero, one}};
  constexpr polynomial_nttp<GF5, 3> a{{four, three, one, two}};
  constexpr auto r = poly_rem<GF5, 2, mod>(a);
  constexpr polynomial_nttp<GF5, 1> q{{one, two}};
  constexpr auto reconstructed = (q * mod) + r;
  return reconstructed[0] == a[0] && reconstructed[1] == a[1] && reconstructed[2] == a[2] && reconstructed[3] == a[3];
}

// ============================================================
// Runtime coverage — volatile-sourced dividend coefficients
// force runtime execution of the same code paths
// ============================================================

bool runtime_test_poly_rem_gf7_cubic_quadratic()
{
  volatile std::uint64_t v3 = 3;
  volatile std::uint64_t v2 = 2;
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 3> a{{GF7(v3), GF7(v2), GF7(v0), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(1)}};
  auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == GF7(3) && r[1] == GF7(1);
}

bool runtime_test_poly_rem_gf7_exact()
{
  volatile std::uint64_t v6 = 6;
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 2> a{{GF7(v6), GF7(v0), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 1> mod{{GF7(1), GF7(1)}};
  auto r = poly_rem<GF7, 1, mod>(a);
  return r[0] == GF7(0);
}

bool runtime_test_poly_rem_gf5_cubic()
{
  volatile std::uint64_t v4 = 4;
  volatile std::uint64_t v3 = 3;
  volatile std::uint64_t v1 = 1;
  volatile std::uint64_t v2 = 2;
  polynomial_nttp<GF5, 3> a{{GF5(v4), GF5(v3), GF5(v1), GF5(v2)}};
  constexpr polynomial_nttp<GF5, 2> mod{{GF5(2), GF5(0), GF5(1)}};
  auto r = poly_rem<GF5, 2, mod>(a);
  return r[0] == GF5(2) && r[1] == GF5(4);
}

bool runtime_test_poly_rem_gf2_cubic()
{
  // x³ mod (x²+x+1) = 1 at runtime
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF2, 3> a{{GF2(v0), GF2(v0), GF2(v0), GF2(v1)}};
  constexpr polynomial_nttp<GF2, 2> mod{{GF2(1), GF2(1), GF2(1)}};
  auto r = poly_rem<GF2, 2, mod>(a);
  return r[0] == GF2(1) && r[1] == GF2(0);
}

bool runtime_test_poly_rem_gf7_nonmonic()
{
  // (x³+x+1) mod (3x²+1) = 3x+1 → {1,3} at runtime
  volatile std::uint64_t v1 = 1;
  volatile std::uint64_t v0 = 0;
  polynomial_nttp<GF7, 3> a{{GF7(v1), GF7(v1), GF7(v0), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(3)}};
  auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == GF7(1) && r[1] == GF7(3);
}

bool runtime_test_poly_rem_gf17_cubic()
{
  // (x³+x+1) mod (x²+3) = 15x+1 → {1,15} at runtime
  volatile std::uint64_t v1 = 1;
  volatile std::uint64_t v0 = 0;
  polynomial_nttp<GF17, 3> a{{GF17(v1), GF17(v1), GF17(v0), GF17(v1)}};
  constexpr polynomial_nttp<GF17, 2> mod{{GF17(3), GF17(0), GF17(1)}};
  auto r = poly_rem<GF17, 2, mod>(a);
  return r[0] == GF17(1) && r[1] == GF17(15);
}

bool runtime_test_poly_rem_gf7_degree4_quintic()
{
  // x⁵+2x+3 mod (x⁴+1) = x+3 → {3,1,0,0} at runtime
  volatile std::uint64_t v3 = 3;
  volatile std::uint64_t v2 = 2;
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 5> a{{GF7(v3), GF7(v2), GF7(v0), GF7(v0), GF7(v0), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 4> mod{{GF7(1), GF7(0), GF7(0), GF7(0), GF7(1)}};
  auto r = poly_rem<GF7, 4, mod>(a);
  return r[0] == GF7(3) && r[1] == GF7(1) && r[2] == GF7(0) && r[3] == GF7(0);
}

bool runtime_test_poly_rem_gf7_zero_dividend()
{
  volatile std::uint64_t v0 = 0;
  polynomial_nttp<GF7, 3> a{{GF7(v0), GF7(v0), GF7(v0), GF7(v0)}};
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(1)}};
  auto r = poly_rem<GF7, 2, mod>(a);
  return r[0] == GF7(0) && r[1] == GF7(0);
}

// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

static_assert(test_poly_rem_gf7_m_lt_n_constant(), "poly_rem GF(7): M=0 < N=2, constant");
static_assert(test_poly_rem_gf7_m_lt_n_linear(), "poly_rem GF(7): M=1 < N=2, linear passthrough");
static_assert(test_poly_rem_gf7_exact_linear(), "poly_rem GF(7): x²-1 mod (x+1) = 0");
static_assert(test_poly_rem_gf7_quadratic_reduces(), "poly_rem GF(7): x²+x+1 mod (x²+1) = x");
static_assert(test_poly_rem_gf7_cubic_quadratic(), "poly_rem GF(7): x³+2x+3 mod (x²+1) = x+3");
static_assert(test_poly_rem_gf7_quartic_quadratic(), "poly_rem GF(7): x⁴+x²+1 mod (x²+1) = 1");
static_assert(test_poly_rem_gf7_quartic_linear_remainder(), "poly_rem GF(7): x⁴+x mod (x²+1) = x+1");
static_assert(test_poly_rem_gf7_zero_dividend(), "poly_rem GF(7): zero dividend");
static_assert(test_poly_rem_gf5_constant_passthrough(), "poly_rem GF(5): constant passthrough");
static_assert(test_poly_rem_gf5_quadratic_reduces(), "poly_rem GF(5): x²+x+1 mod (x²+2) = x+4");
static_assert(test_poly_rem_gf5_cubic_quadratic(), "poly_rem GF(5): 2x³+x²+3x+4 mod (x²+2) = 4x+2");
static_assert(test_poly_rem_gf7_cubic_mod_m_lt_n(), "poly_rem GF(7): deg-1 mod deg-3, passthrough");
static_assert(test_poly_rem_gf7_cubic_mod_degree3(), "poly_rem GF(7): x³+1 mod (x³+x+1) = 6x");
static_assert(test_poly_rem_gf2_constant_passthrough(), "poly_rem GF(2): constant passthrough");
static_assert(test_poly_rem_gf2_linear_passthrough(), "poly_rem GF(2): linear passthrough");
static_assert(test_poly_rem_gf2_x2_plus_x(), "poly_rem GF(2): x²+x mod (x²+x+1) = 1");
static_assert(test_poly_rem_gf2_cubic(), "poly_rem GF(2): x³ mod (x²+x+1) = 1");
static_assert(test_poly_rem_gf2_quartic(), "poly_rem GF(2): x⁴ mod (x²+x+1) = x");
static_assert(test_poly_rem_gf2_reconstruction(), "poly_rem GF(2): reconstruction a = q*mod + r");
static_assert(test_poly_rem_gf7_nonmonic_cubic(), "poly_rem GF(7): nonmonic 3x²+1, cubic dividend");
static_assert(test_poly_rem_gf7_nonmonic_m_equals_n(), "poly_rem GF(7): nonmonic 3x²+1, M=N");
static_assert(test_poly_rem_gf7_nonmonic_exact(), "poly_rem GF(7): nonmonic 3x²+1, exact divisibility");
static_assert(test_poly_rem_gf17_constant_passthrough(), "poly_rem GF(17): constant passthrough");
static_assert(test_poly_rem_gf17_cubic(), "poly_rem GF(17): x³+x+1 mod (x²+3) = 15x+1");
static_assert(test_poly_rem_gf17_quartic(), "poly_rem GF(17): x⁴+2x+5 mod (x²+3) = 2x+14");
static_assert(test_poly_rem_gf17_reconstruction(), "poly_rem GF(17): reconstruction a = q*mod + r");
static_assert(test_poly_rem_gf7_degree4_mod_quintic(), "poly_rem GF(7): deg-4 mod, x⁵+2x+3 → x+3");
static_assert(test_poly_rem_gf7_degree4_mod_x7(), "poly_rem GF(7): deg-4 mod, x⁷+1 → 6x³+1");
static_assert(test_poly_rem_gf7_degree4_mod_exact_multiple(), "poly_rem GF(7): deg-4 mod, exact multiple");
static_assert(test_poly_rem_gf5_reconstruction(), "poly_rem GF(5): reconstruction a = q*mod + r");

int main()
{
  constexpr bool ct_gf7_const = test_poly_rem_gf7_m_lt_n_constant();
  constexpr bool ct_gf7_linear = test_poly_rem_gf7_m_lt_n_linear();
  constexpr bool ct_gf7_exact = test_poly_rem_gf7_exact_linear();
  constexpr bool ct_gf7_quad = test_poly_rem_gf7_quadratic_reduces();
  constexpr bool ct_gf7_cubic = test_poly_rem_gf7_cubic_quadratic();
  constexpr bool ct_gf7_quartic = test_poly_rem_gf7_quartic_quadratic();
  constexpr bool ct_gf7_qlinrem = test_poly_rem_gf7_quartic_linear_remainder();
  constexpr bool ct_gf7_zero = test_poly_rem_gf7_zero_dividend();
  constexpr bool ct_gf5_const = test_poly_rem_gf5_constant_passthrough();
  constexpr bool ct_gf5_quad = test_poly_rem_gf5_quadratic_reduces();
  constexpr bool ct_gf5_cubic = test_poly_rem_gf5_cubic_quadratic();
  constexpr bool ct_gf7_cm_lt_n = test_poly_rem_gf7_cubic_mod_m_lt_n();
  constexpr bool ct_gf7_cd3 = test_poly_rem_gf7_cubic_mod_degree3();
  constexpr bool ct_gf2_const = test_poly_rem_gf2_constant_passthrough();
  constexpr bool ct_gf2_linear = test_poly_rem_gf2_linear_passthrough();
  constexpr bool ct_gf2_x2x = test_poly_rem_gf2_x2_plus_x();
  constexpr bool ct_gf2_cubic = test_poly_rem_gf2_cubic();
  constexpr bool ct_gf2_quartic = test_poly_rem_gf2_quartic();
  constexpr bool ct_gf2_recon = test_poly_rem_gf2_reconstruction();
  constexpr bool ct_nm_cubic = test_poly_rem_gf7_nonmonic_cubic();
  constexpr bool ct_nm_mn = test_poly_rem_gf7_nonmonic_m_equals_n();
  constexpr bool ct_nm_exact = test_poly_rem_gf7_nonmonic_exact();
  constexpr bool ct_gf17_const = test_poly_rem_gf17_constant_passthrough();
  constexpr bool ct_gf17_cubic = test_poly_rem_gf17_cubic();
  constexpr bool ct_gf17_quartic = test_poly_rem_gf17_quartic();
  constexpr bool ct_gf17_recon = test_poly_rem_gf17_reconstruction();
  constexpr bool ct_d4_quintic = test_poly_rem_gf7_degree4_mod_quintic();
  constexpr bool ct_d4_x7 = test_poly_rem_gf7_degree4_mod_x7();
  constexpr bool ct_d4_exact = test_poly_rem_gf7_degree4_mod_exact_multiple();
  constexpr bool ct_gf7_recon = test_poly_rem_gf7_reconstruction();
  constexpr bool ct_gf5_recon = test_poly_rem_gf5_reconstruction();

  bool rt_gf7_cubic = runtime_test_poly_rem_gf7_cubic_quadratic();
  bool rt_gf7_exact = runtime_test_poly_rem_gf7_exact();
  bool rt_gf5_cubic = runtime_test_poly_rem_gf5_cubic();
  bool rt_gf2_cubic = runtime_test_poly_rem_gf2_cubic();
  bool rt_nm = runtime_test_poly_rem_gf7_nonmonic();
  bool rt_gf17_cubic = runtime_test_poly_rem_gf17_cubic();
  bool rt_d4_quintic = runtime_test_poly_rem_gf7_degree4_quintic();
  bool rt_zero = runtime_test_poly_rem_gf7_zero_dividend();

  if constexpr (ct_gf7_const && ct_gf7_linear && ct_gf7_exact && ct_gf7_quad && ct_gf7_cubic && ct_gf7_quartic &&
                ct_gf7_qlinrem && ct_gf7_zero && ct_gf5_const && ct_gf5_quad && ct_gf5_cubic && ct_gf7_cm_lt_n &&
                ct_gf7_cd3 && ct_gf2_const && ct_gf2_linear && ct_gf2_x2x && ct_gf2_cubic && ct_gf2_quartic &&
                ct_gf2_recon && ct_nm_cubic && ct_nm_mn && ct_nm_exact && ct_gf17_const && ct_gf17_cubic &&
                ct_gf17_quartic && ct_gf17_recon && ct_d4_quintic && ct_d4_x7 && ct_d4_exact && ct_gf7_recon &&
                ct_gf5_recon)
  {
    if (rt_gf7_cubic && rt_gf7_exact && rt_gf5_cubic && rt_gf2_cubic && rt_nm && rt_gf17_cubic && rt_d4_quintic &&
        rt_zero)
      return 0;
    return 1;
  }
  else
  {
    return 1;
  }
}
