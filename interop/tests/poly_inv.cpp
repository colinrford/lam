/*
 *  poly_inv.cpp – Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit tests for poly_inv: multiplicative inverse in R[x]/(Mod)
 *
 *  Tests are structured as consteval functions verified by static_assert
 *  (compile-time) and also exercised with volatile-sourced input coefficients
 *  to force runtime execution (runtime coverage).
 *
 *  Tests resistant to compile-time evaluation: none; all ZqElement arithmetic
 *    (including / which calls mod_inv) is constexpr. The XGCD loop terminates
 *    within a bounded number of steps for small degrees, but the cumulative
 *    constexpr step count can be large — compile with -fconstexpr-steps=100000000.
 *  Tests resistant to runtime: none — Mod is a NTTP (compile-time), but input
 *    polynomial coefficients can be volatile-sourced.
 *
 *  Quotient rings exercised:
 *    GF(4)   = GF(2)[x]/(x²+x+1)   — characteristic 2; -x = x, only two non-trivial elements
 *    GF(49)  = GF(7)[x]/(x²+1)      — odd prime, degree-2 irreducible
 *    GF(25)  = GF(5)[x]/(x²+2)      — different odd prime
 *    GF(289) = GF(17)[x]/(x²+3)     — larger prime; x²+3 irreducible over GF(17)
 *    GF(7³)  = GF(7)[x]/(x³+x+1)    — degree-3, exercises more XGCD iterations
 *
 *  All expected inverse values are hand-verified by solving the linear system
 *  arising from polynomial multiplication mod the irreducible.
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
// File-scope moduli (constexpr, usable as NTTPs throughout)
// ============================================================

constexpr polynomial_nttp<GF2, 2> gf4_mod{{GF2(1), GF2(1), GF2(1)}};           // x²+x+1 over GF(2)
constexpr polynomial_nttp<GF7, 2> gf49_mod{{GF7(1), GF7(0), GF7(1)}};          // x²+1 over GF(7)
constexpr polynomial_nttp<GF5, 2> gf25_mod{{GF5(2), GF5(0), GF5(1)}};          // x²+2 over GF(5)
constexpr polynomial_nttp<GF17, 2> gf289_mod{{GF17(3), GF17(0), GF17(1)}};     // x²+3 over GF(17)
constexpr polynomial_nttp<GF7, 3> gf7_3_mod{{GF7(1), GF7(1), GF7(0), GF7(1)}}; // x³+x+1

// ============================================================
// GF(4) = GF(2)[x]/(x²+x+1)
//
// Only two non-trivial, non-unity elements: x and x+1.
// In characteristic 2: -x = x, 1+1 = 0, so subtraction = addition.
//
// Inverse derivation (solve (cx+d)(ax+b) ≡ 1 mod x²+x+1):
//   (ax+b)(cx+d) = acx²+(ad+bc)x+bd
//   x² ≡ x+1 (from x²+x+1=0 over GF(2), so x²=-x-1=x+1)
//   ≡ ac(x+1)+(ad+bc)x+bd = (ac+ad+bc)x+(ac+bd)
//
//   inv(x):   c=1,d=0 → (a+b)x+a = 1; a=1, a+b=0 → b=1  → inv(x) = x+1 = {1,1}
//   inv(x+1): c=1,d=1 → (a+a+b)x+(a+b) = bx+(a+b) = 1; b=0, a=1 → inv(x+1) = x = {0,1}
//
// Note: x and x+1 are inverses of each other, which is verified by roundtrip.
// ============================================================

consteval bool test_poly_inv_gf4_x()
{
  // inv(x) in GF(4) = x+1 = {1,1}
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 1> a{{zero, one}}; // x
  constexpr auto inv_a = poly_inv<GF2, 2, gf4_mod>(a);
  return inv_a[0] == one && inv_a[1] == one;
}

consteval bool test_poly_inv_gf4_x_plus_1()
{
  // inv(x+1) in GF(4) = x = {0,1}
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 1> a{{one, one}}; // x+1
  constexpr auto inv_a = poly_inv<GF2, 2, gf4_mod>(a);
  return inv_a[0] == zero && inv_a[1] == one;
}

consteval bool test_poly_inv_gf4_roundtrip_x()
{
  // x * (x+1) = x²+x ≡ (x+1)+x = 1 in GF(2) (x²≡x+1)
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 1> a{{zero, one}};    // x
  constexpr polynomial_nttp<GF2, 1> inv_a{{one, one}}; // x+1
  constexpr auto product = a * inv_a;
  constexpr auto r = poly_rem<GF2, 2, gf4_mod>(product);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_inv_gf4_roundtrip_self_consistent_x()
{
  // Compute inv(x) via poly_inv and verify a·inv(a) ≡ 1
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 1> a{{zero, one}};
  constexpr auto inv_a = poly_inv<GF2, 2, gf4_mod>(a);
  constexpr auto r = poly_rem<GF2, 2, gf4_mod>(a * inv_a);
  return r[0] == one && r[1] == zero;
}

consteval bool test_poly_inv_gf4_roundtrip_self_consistent_x_plus_1()
{
  constexpr GF2 zero(0);
  constexpr GF2 one(1);
  constexpr polynomial_nttp<GF2, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF2, 2, gf4_mod>(a);
  constexpr auto r = poly_rem<GF2, 2, gf4_mod>(a * inv_a);
  return r[0] == one && r[1] == zero;
}

// ============================================================
// GF(49) = GF(7)[x]/(x²+1)
//
// Inverse derivation (x² ≡ -1):
//   (ax+b)·(cx+d) ≡ -ac+(ad+bc)x+bd
//   inv(x):    c=1,d=0 → -a=1,b=0 → a=6              → {0,6}
//   inv(x+1):  c=1,d=1 → -a+b=1, a+b=0 → b=-a, -2a=1, a=3, b=4 → {4,3}
//   inv(x+2):  c=1,d=2 → -a+2b=1, 2a+b=0 → b=-2a, -5a=1, a=4, b=6 → {6,4}
//   inv(3x+5): c=3,d=5 → -3a+5b=1, 5a+3b=0 → b=-5a/3=3a; -3a+15a=1, 12a=1, 5a=1, a=3, b=2 → {2,3}
// ============================================================

consteval bool test_poly_inv_gf49_constant_one()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 1> a{{one, zero}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  return inv_a[0] == one && inv_a[1] == zero;
}

consteval bool test_poly_inv_gf49_x()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 1> a{{zero, one}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  return inv_a[0] == zero && inv_a[1] == six;
}

consteval bool test_poly_inv_gf49_x_plus_1()
{
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr GF7 four(4);
  constexpr polynomial_nttp<GF7, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  return inv_a[0] == four && inv_a[1] == three;
}

consteval bool test_poly_inv_gf49_x_plus_2()
{
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 four(4);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 1> a{{two, one}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  return inv_a[0] == six && inv_a[1] == four;
}

consteval bool test_poly_inv_gf49_3x_plus_5()
{
  constexpr GF7 two(2);
  constexpr GF7 three(3);
  constexpr GF7 five(5);
  constexpr polynomial_nttp<GF7, 1> a{{five, three}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  return inv_a[0] == two && inv_a[1] == three;
}

consteval bool test_poly_inv_gf49_roundtrip_x()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 1> a{{zero, one}};
  constexpr polynomial_nttp<GF7, 1> inv_a{{zero, six}};
  constexpr auto r = poly_rem<GF7, 2, gf49_mod>(a * inv_a);
  return r[0] == one && r[1] == GF7(0);
}

consteval bool test_poly_inv_gf49_roundtrip_x_plus_1()
{
  constexpr GF7 one(1);
  constexpr GF7 three(3);
  constexpr GF7 four(4);
  constexpr polynomial_nttp<GF7, 1> a{{one, one}};
  constexpr polynomial_nttp<GF7, 1> inv_a{{four, three}};
  constexpr auto r = poly_rem<GF7, 2, gf49_mod>(a * inv_a);
  return r[0] == one && r[1] == GF7(0);
}

consteval bool test_poly_inv_gf49_roundtrip_x_plus_2()
{
  constexpr GF7 one(1);
  constexpr GF7 two(2);
  constexpr GF7 four(4);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 1> a{{two, one}};
  constexpr polynomial_nttp<GF7, 1> inv_a{{six, four}};
  constexpr auto r = poly_rem<GF7, 2, gf49_mod>(a * inv_a);
  return r[0] == one && r[1] == GF7(0);
}

consteval bool test_poly_inv_gf49_roundtrip_self_consistent()
{
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF7, 2, gf49_mod>(a);
  constexpr auto r = poly_rem<GF7, 2, gf49_mod>(a * inv_a);
  return r[0] == one && r[1] == GF7(0);
}

// ============================================================
// GF(25) = GF(5)[x]/(x²+2)
//
// x² ≡ -2 = 3 in GF(5):
//   inv(x):     -2a=1, b=0 → a=2                 → {0,2}
//   inv(x+1):   -2a+b=1, a+b=0 → b=-a, -3a=1, a=3, b=2 → {2,3}
//   inv(2x+3):  -4a+3b=1, 3a+2b=0 → b=-3a/2=a; -4a+3a=1, -a=1, a=4, b=4 → {4,4}
// ============================================================

consteval bool test_poly_inv_gf25_x()
{
  constexpr GF5 zero(0);
  constexpr GF5 one(1);
  constexpr GF5 two(2);
  constexpr polynomial_nttp<GF5, 1> a{{zero, one}};
  constexpr auto inv_a = poly_inv<GF5, 2, gf25_mod>(a);
  return inv_a[0] == zero && inv_a[1] == two;
}

consteval bool test_poly_inv_gf25_x_plus_1()
{
  constexpr GF5 one(1);
  constexpr GF5 two(2);
  constexpr GF5 three(3);
  constexpr polynomial_nttp<GF5, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF5, 2, gf25_mod>(a);
  return inv_a[0] == two && inv_a[1] == three;
}

consteval bool test_poly_inv_gf25_2x_plus_3()
{
  constexpr GF5 two(2);
  constexpr GF5 three(3);
  constexpr GF5 four(4);
  constexpr polynomial_nttp<GF5, 1> a{{three, two}};
  constexpr auto inv_a = poly_inv<GF5, 2, gf25_mod>(a);
  return inv_a[0] == four && inv_a[1] == four;
}

consteval bool test_poly_inv_gf25_roundtrip_x()
{
  constexpr GF5 zero(0);
  constexpr GF5 one(1);
  constexpr GF5 two(2);
  constexpr polynomial_nttp<GF5, 1> a{{zero, one}};
  constexpr polynomial_nttp<GF5, 1> inv_a{{zero, two}};
  constexpr auto r = poly_rem<GF5, 2, gf25_mod>(a * inv_a);
  return r[0] == one && r[1] == GF5(0);
}

consteval bool test_poly_inv_gf25_roundtrip_self_consistent()
{
  constexpr GF5 one(1);
  constexpr polynomial_nttp<GF5, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF5, 2, gf25_mod>(a);
  constexpr auto r = poly_rem<GF5, 2, gf25_mod>(a * inv_a);
  return r[0] == one && r[1] == GF5(0);
}

// ============================================================
// GF(289) = GF(17)[x]/(x²+3)
//
// x² ≡ -3 = 14 in GF(17):
//   inv(x):    -3a=1, b=0 → a=-1/3; 3⁻¹=6, a=-6=11     → {0,11}
//   inv(x+1):  -3a+b=1, a+b=0 → b=-a, -4a=1; 4⁻¹=13, a=-13=4, b=13 → {13,4}
//   inv(3x+5): -9a+5b=1, 5a+3b=0 → b=-5a/3=-10a=7a; -9a+35a=1, 26a=1=9a (26≡9),
//              9⁻¹ mod 17: 9·2=18≡1, so 9⁻¹=2; a=2, b=7·2=14; → {14,6}
//              Wait, let me redo this: 5a+3b=0 → b=-5a/3; 3⁻¹=6, b=-5a·6=-30a=-13a=4a.
//              -9a+5(4a)=1 → -9a+20a=1 → 11a=1; 11⁻¹ mod 17: 11·14=154=9·17+1, so 11⁻¹=14.
//              a=14, b=4·14=56=56-3·17=56-51=5 → {5,14}
// ============================================================

consteval bool test_poly_inv_gf289_x()
{
  // inv(x) in GF(289): 11x → {0,11}
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 eleven(11);
  constexpr polynomial_nttp<GF17, 1> a{{zero, one}};
  constexpr auto inv_a = poly_inv<GF17, 2, gf289_mod>(a);
  return inv_a[0] == zero && inv_a[1] == eleven;
}

consteval bool test_poly_inv_gf289_x_plus_1()
{
  // inv(x+1) in GF(289): 4x+13 → {13,4}
  constexpr GF17 one(1);
  constexpr GF17 four(4);
  constexpr GF17 thirteen(13);
  constexpr polynomial_nttp<GF17, 1> a{{one, one}};
  constexpr auto inv_a = poly_inv<GF17, 2, gf289_mod>(a);
  return inv_a[0] == thirteen && inv_a[1] == four;
}

consteval bool test_poly_inv_gf289_3x_plus_5()
{
  // inv(3x+5) in GF(289): 14x+5 → {5,14}
  constexpr GF17 five(5);
  constexpr GF17 fourteen(14);
  constexpr polynomial_nttp<GF17, 1> a{{five, GF17(3)}};
  constexpr auto inv_a = poly_inv<GF17, 2, gf289_mod>(a);
  return inv_a[0] == five && inv_a[1] == fourteen;
}

consteval bool test_poly_inv_gf289_roundtrip_x()
{
  constexpr GF17 zero(0);
  constexpr GF17 one(1);
  constexpr GF17 eleven(11);
  constexpr polynomial_nttp<GF17, 1> a{{zero, one}};
  constexpr polynomial_nttp<GF17, 1> inv_a{{zero, eleven}};
  constexpr auto r = poly_rem<GF17, 2, gf289_mod>(a * inv_a);
  return r[0] == one && r[1] == GF17(0);
}

consteval bool test_poly_inv_gf289_roundtrip_x_plus_1()
{
  constexpr GF17 one(1);
  constexpr GF17 four(4);
  constexpr GF17 thirteen(13);
  constexpr polynomial_nttp<GF17, 1> a{{one, one}};
  constexpr polynomial_nttp<GF17, 1> inv_a{{thirteen, four}};
  constexpr auto r = poly_rem<GF17, 2, gf289_mod>(a * inv_a);
  return r[0] == one && r[1] == GF17(0);
}

consteval bool test_poly_inv_gf289_roundtrip_self_consistent()
{
  constexpr GF17 five(5);
  constexpr polynomial_nttp<GF17, 1> a{{five, GF17(3)}}; // 3x+5
  constexpr auto inv_a = poly_inv<GF17, 2, gf289_mod>(a);
  constexpr auto r = poly_rem<GF17, 2, gf289_mod>(a * inv_a);
  return r[0] == GF17(1) && r[1] == GF17(0);
}

// ============================================================
// GF(7³) = GF(7)[x]/(x³+x+1)
//
// Inverse derivation for a = x (x³ ≡ -x-1):
//   x·(ax²+bx+c) = ax³+bx²+cx ≡ -a(x+1)+bx²+cx = bx²+(c-a)x+(-a)
//   b=0, c-a=0, -a=1 → a=6, c=6 → inv(x) = 6x²+6 = {6,0,6}
// ============================================================

consteval bool test_poly_inv_gf7_3_x()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 2> a{{zero, one, zero}};
  constexpr auto inv_a = poly_inv<GF7, 3, gf7_3_mod>(a);
  return inv_a[0] == six && inv_a[1] == zero && inv_a[2] == six;
}

consteval bool test_poly_inv_gf7_3_roundtrip_x()
{
  constexpr GF7 zero(0);
  constexpr GF7 one(1);
  constexpr GF7 six(6);
  constexpr polynomial_nttp<GF7, 2> a{{zero, one, zero}};
  constexpr polynomial_nttp<GF7, 2> inv_a{{six, zero, six}};
  constexpr auto r = poly_rem<GF7, 3, gf7_3_mod>(a * inv_a);
  return r[0] == one && r[1] == zero && r[2] == zero;
}

consteval bool test_poly_inv_gf7_3_roundtrip_self_consistent()
{
  constexpr GF7 one(1);
  constexpr polynomial_nttp<GF7, 2> a{{one, one, one}}; // x²+x+1
  constexpr auto inv_a = poly_inv<GF7, 3, gf7_3_mod>(a);
  constexpr auto r = poly_rem<GF7, 3, gf7_3_mod>(a * inv_a);
  return r[0] == one && r[1] == GF7(0) && r[2] == GF7(0);
}

// ============================================================
// Runtime coverage — volatile-sourced input polynomial coefficients
// ============================================================

bool runtime_test_poly_inv_gf4_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF2, 1> a{{GF2(v0), GF2(v1)}};
  constexpr polynomial_nttp<GF2, 2> mod{{GF2(1), GF2(1), GF2(1)}};
  auto inv_a = poly_inv<GF2, 2, mod>(a);
  return inv_a[0] == GF2(1) && inv_a[1] == GF2(1);
}

bool runtime_test_poly_inv_gf4_x_plus_1()
{
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF2, 1> a{{GF2(v1), GF2(v1)}};
  constexpr polynomial_nttp<GF2, 2> mod{{GF2(1), GF2(1), GF2(1)}};
  auto inv_a = poly_inv<GF2, 2, mod>(a);
  return inv_a[0] == GF2(0) && inv_a[1] == GF2(1);
}

bool runtime_test_poly_inv_gf49_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 1> a{{GF7(v0), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(1)}};
  auto inv_a = poly_inv<GF7, 2, mod>(a);
  return inv_a[0] == GF7(0) && inv_a[1] == GF7(6);
}

bool runtime_test_poly_inv_gf49_x_plus_1()
{
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 1> a{{GF7(v1), GF7(v1)}};
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(1)}};
  auto inv_a = poly_inv<GF7, 2, mod>(a);
  return inv_a[0] == GF7(4) && inv_a[1] == GF7(3);
}

bool runtime_test_poly_inv_gf49_roundtrip()
{
  // Compute inv at runtime and verify the product is 1 mod gf49_mod
  volatile std::uint64_t v2 = 2;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 1> a{{GF7(v2), GF7(v1)}}; // x + 2
  constexpr polynomial_nttp<GF7, 2> mod{{GF7(1), GF7(0), GF7(1)}};
  auto inv_a = poly_inv<GF7, 2, mod>(a);
  auto r = poly_rem<GF7, 2, mod>(a * inv_a);
  return r[0] == GF7(1) && r[1] == GF7(0);
}

bool runtime_test_poly_inv_gf25_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF5, 1> a{{GF5(v0), GF5(v1)}};
  constexpr polynomial_nttp<GF5, 2> mod{{GF5(2), GF5(0), GF5(1)}};
  auto inv_a = poly_inv<GF5, 2, mod>(a);
  return inv_a[0] == GF5(0) && inv_a[1] == GF5(2);
}

bool runtime_test_poly_inv_gf289_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF17, 1> a{{GF17(v0), GF17(v1)}};
  constexpr polynomial_nttp<GF17, 2> mod{{GF17(3), GF17(0), GF17(1)}};
  auto inv_a = poly_inv<GF17, 2, mod>(a);
  return inv_a[0] == GF17(0) && inv_a[1] == GF17(11);
}

bool runtime_test_poly_inv_gf7_3_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  polynomial_nttp<GF7, 2> a{{GF7(v0), GF7(v1), GF7(v0)}};
  constexpr polynomial_nttp<GF7, 3> mod{{GF7(1), GF7(1), GF7(0), GF7(1)}};
  auto inv_a = poly_inv<GF7, 3, mod>(a);
  return inv_a[0] == GF7(6) && inv_a[1] == GF7(0) && inv_a[2] == GF7(6);
}

// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

static_assert(test_poly_inv_gf4_x(), "poly_inv GF(4): inv(x) = x+1");
static_assert(test_poly_inv_gf4_x_plus_1(), "poly_inv GF(4): inv(x+1) = x");
static_assert(test_poly_inv_gf4_roundtrip_x(), "poly_inv GF(4): x*(x+1) = 1");
static_assert(test_poly_inv_gf4_roundtrip_self_consistent_x(), "poly_inv GF(4): self-consistent roundtrip, x");
static_assert(test_poly_inv_gf4_roundtrip_self_consistent_x_plus_1(), "poly_inv GF(4): self-consistent roundtrip, x+1");
static_assert(test_poly_inv_gf49_constant_one(), "poly_inv GF(49): inv(1) = 1");
static_assert(test_poly_inv_gf49_x(), "poly_inv GF(49): inv(x) = 6x");
static_assert(test_poly_inv_gf49_x_plus_1(), "poly_inv GF(49): inv(x+1) = 3x+4");
static_assert(test_poly_inv_gf49_x_plus_2(), "poly_inv GF(49): inv(x+2) = 4x+6");
static_assert(test_poly_inv_gf49_3x_plus_5(), "poly_inv GF(49): inv(3x+5) = 3x+2");
static_assert(test_poly_inv_gf49_roundtrip_x(), "poly_inv GF(49): x*inv(x) = 1");
static_assert(test_poly_inv_gf49_roundtrip_x_plus_1(), "poly_inv GF(49): (x+1)*inv(x+1) = 1");
static_assert(test_poly_inv_gf49_roundtrip_x_plus_2(), "poly_inv GF(49): (x+2)*inv(x+2) = 1");
static_assert(test_poly_inv_gf49_roundtrip_self_consistent(), "poly_inv GF(49): self-consistent roundtrip");
static_assert(test_poly_inv_gf25_x(), "poly_inv GF(25): inv(x) = 2x");
static_assert(test_poly_inv_gf25_x_plus_1(), "poly_inv GF(25): inv(x+1) = 3x+2");
static_assert(test_poly_inv_gf25_2x_plus_3(), "poly_inv GF(25): inv(2x+3) = 4x+4");
static_assert(test_poly_inv_gf25_roundtrip_x(), "poly_inv GF(25): x*inv(x) = 1");
static_assert(test_poly_inv_gf25_roundtrip_self_consistent(), "poly_inv GF(25): self-consistent roundtrip");
static_assert(test_poly_inv_gf289_x(), "poly_inv GF(289): inv(x) = 11x");
static_assert(test_poly_inv_gf289_x_plus_1(), "poly_inv GF(289): inv(x+1) = 4x+13");
static_assert(test_poly_inv_gf289_3x_plus_5(), "poly_inv GF(289): inv(3x+5) = 14x+5");
static_assert(test_poly_inv_gf289_roundtrip_x(), "poly_inv GF(289): x*inv(x) = 1");
static_assert(test_poly_inv_gf289_roundtrip_x_plus_1(), "poly_inv GF(289): (x+1)*inv(x+1) = 1");
static_assert(test_poly_inv_gf289_roundtrip_self_consistent(), "poly_inv GF(289): self-consistent roundtrip");
static_assert(test_poly_inv_gf7_3_x(), "poly_inv GF(7³): inv(x) = 6x²+6");
static_assert(test_poly_inv_gf7_3_roundtrip_x(), "poly_inv GF(7³): x*inv(x) = 1");
static_assert(test_poly_inv_gf7_3_roundtrip_self_consistent(), "poly_inv GF(7³): self-consistent roundtrip");

int main()
{
  constexpr bool ct_gf4_x = test_poly_inv_gf4_x();
  constexpr bool ct_gf4_xp1 = test_poly_inv_gf4_x_plus_1();
  constexpr bool ct_gf4_rt_x = test_poly_inv_gf4_roundtrip_x();
  constexpr bool ct_gf4_sc_x = test_poly_inv_gf4_roundtrip_self_consistent_x();
  constexpr bool ct_gf4_sc_xp1 = test_poly_inv_gf4_roundtrip_self_consistent_x_plus_1();
  constexpr bool ct_gf49_one = test_poly_inv_gf49_constant_one();
  constexpr bool ct_gf49_x = test_poly_inv_gf49_x();
  constexpr bool ct_gf49_xp1 = test_poly_inv_gf49_x_plus_1();
  constexpr bool ct_gf49_xp2 = test_poly_inv_gf49_x_plus_2();
  constexpr bool ct_gf49_3xp5 = test_poly_inv_gf49_3x_plus_5();
  constexpr bool ct_gf49_rt_x = test_poly_inv_gf49_roundtrip_x();
  constexpr bool ct_gf49_rt_xp1 = test_poly_inv_gf49_roundtrip_x_plus_1();
  constexpr bool ct_gf49_rt_xp2 = test_poly_inv_gf49_roundtrip_x_plus_2();
  constexpr bool ct_gf49_rt_sc = test_poly_inv_gf49_roundtrip_self_consistent();
  constexpr bool ct_gf25_x = test_poly_inv_gf25_x();
  constexpr bool ct_gf25_xp1 = test_poly_inv_gf25_x_plus_1();
  constexpr bool ct_gf25_2xp3 = test_poly_inv_gf25_2x_plus_3();
  constexpr bool ct_gf25_rt_x = test_poly_inv_gf25_roundtrip_x();
  constexpr bool ct_gf25_rt_sc = test_poly_inv_gf25_roundtrip_self_consistent();
  constexpr bool ct_gf289_x = test_poly_inv_gf289_x();
  constexpr bool ct_gf289_xp1 = test_poly_inv_gf289_x_plus_1();
  constexpr bool ct_gf289_3xp5 = test_poly_inv_gf289_3x_plus_5();
  constexpr bool ct_gf289_rt_x = test_poly_inv_gf289_roundtrip_x();
  constexpr bool ct_gf289_rt_xp1 = test_poly_inv_gf289_roundtrip_x_plus_1();
  constexpr bool ct_gf289_rt_sc = test_poly_inv_gf289_roundtrip_self_consistent();
  constexpr bool ct_gf7_3_x = test_poly_inv_gf7_3_x();
  constexpr bool ct_gf7_3_rt_x = test_poly_inv_gf7_3_roundtrip_x();
  constexpr bool ct_gf7_3_rt_sc = test_poly_inv_gf7_3_roundtrip_self_consistent();

  bool rt_gf4_x = runtime_test_poly_inv_gf4_x();
  bool rt_gf4_xp1 = runtime_test_poly_inv_gf4_x_plus_1();
  bool rt_gf49_x = runtime_test_poly_inv_gf49_x();
  bool rt_gf49_xp1 = runtime_test_poly_inv_gf49_x_plus_1();
  bool rt_gf49_rt = runtime_test_poly_inv_gf49_roundtrip();
  bool rt_gf25_x = runtime_test_poly_inv_gf25_x();
  bool rt_gf289_x = runtime_test_poly_inv_gf289_x();
  bool rt_gf7_3_x = runtime_test_poly_inv_gf7_3_x();

  if constexpr (ct_gf4_x && ct_gf4_xp1 && ct_gf4_rt_x && ct_gf4_sc_x && ct_gf4_sc_xp1 && ct_gf49_one && ct_gf49_x &&
                ct_gf49_xp1 && ct_gf49_xp2 && ct_gf49_3xp5 && ct_gf49_rt_x && ct_gf49_rt_xp1 && ct_gf49_rt_xp2 &&
                ct_gf49_rt_sc && ct_gf25_x && ct_gf25_xp1 && ct_gf25_2xp3 && ct_gf25_rt_x && ct_gf25_rt_sc &&
                ct_gf289_x && ct_gf289_xp1 && ct_gf289_3xp5 && ct_gf289_rt_x && ct_gf289_rt_xp1 && ct_gf289_rt_sc &&
                ct_gf7_3_x && ct_gf7_3_rt_x && ct_gf7_3_rt_sc)
  {
    if (rt_gf4_x && rt_gf4_xp1 && rt_gf49_x && rt_gf49_xp1 && rt_gf49_rt && rt_gf25_x && rt_gf289_x && rt_gf7_3_x)
      return 0;
    return 1;
  }
  else
  {
    return 1;
  }
}
