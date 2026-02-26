/*
 *  finite_field_extension_test.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit tests for finite_field_extension<K, N, Irreducible>
 *
 *  Tests are structured as consteval functions verified by static_assert
 *  (compile-time) and also exercised with volatile-sourced input to
 *  force runtime execution.
 *
 *  Tests cover the following quotient rings / operations:
 *    GF(49)  = GF(7)[x]/(x²+1)   — standard odd prime; full suite of operations
 *    GF(4)   = GF(2)[x]/(x²+x+1) — characteristic 2: -x = x, x² = x+1
 *    GF(289) = GF(17)[x]/(x²+3)  — larger prime; x² ≡ 14 (mod 17)
 *    GF(25)  = GF(5)[x]/(x²+2)   — different odd prime
 *
 *  Key values (all hand-verified):
 *    GF(49):  x*x={6,0}, inv(x)={0,6}, inv(x+1)={4,3}, x/(x+1)={4,4}
 *    GF(4):   x*x=x+1, inv(x)=x+1, inv(x+1)=x  (char-2)
 *    GF(289): x*x={14,0}, inv(x)={0,11}, inv(x+1)={13,4}, inv(3x+5)={5,14}
 *    GF(25):  x*x={3,0}, inv(x)={0,2}
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.concepts;

using namespace lam;
using namespace lam::cbn;

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)

using GF2 = ZqElement<std::uint64_t, 2>;
using GF5 = ZqElement<std::uint64_t, 5>;
using GF7 = ZqElement<std::uint64_t, 7>;
using GF17 = ZqElement<std::uint64_t, 17>;

// File-scope moduli (constexpr, usable as NTTPs in type aliases below)
constexpr polynomial_nttp<GF2, 2> gf4_mod{{GF2(1), GF2(1), GF2(1)}};       // x²+x+1 over GF(2)
constexpr polynomial_nttp<GF7, 2> gf49_mod{{GF7(1), GF7(0), GF7(1)}};      // x²+1 over GF(7)
constexpr polynomial_nttp<GF5, 2> gf25_mod{{GF5(2), GF5(0), GF5(1)}};      // x²+2 over GF(5)
constexpr polynomial_nttp<GF17, 2> gf289_mod{{GF17(3), GF17(0), GF17(1)}}; // x²+3 over GF(17)

// Element types for each quotient ring
using GF4 = finite_field_extension<GF2, 2, gf4_mod>;
using GF49 = finite_field_extension<GF7, 2, gf49_mod>;
using GF25 = finite_field_extension<GF5, 2, gf25_mod>;
using GF289 = finite_field_extension<GF17, 2, gf289_mod>;

// Verify all four types satisfy field_element_c_weak
static_assert(lam::concepts::experimental::field_element_c_weak<GF4>, "GF4 must satisfy field_element_c_weak");
static_assert(lam::concepts::experimental::field_element_c_weak<GF49>, "GF49 must satisfy field_element_c_weak");
static_assert(lam::concepts::experimental::field_element_c_weak<GF25>, "GF25 must satisfy field_element_c_weak");
static_assert(lam::concepts::experimental::field_element_c_weak<GF289>, "GF289 must satisfy field_element_c_weak");

// ============================================================
// GF(49) = GF(7)[x]/(x²+1)
// x² ≡ -1 ≡ 6 (mod 7)
// ============================================================

consteval bool test_gf49_zero()
{
  constexpr auto z = GF49::zero();
  return z[0] == GF7(0) && z[1] == GF7(0);
}

consteval bool test_gf49_one()
{
  constexpr auto o = GF49::one();
  return o[0] == GF7(1) && o[1] == GF7(0);
}

consteval bool test_gf49_add()
{
  // x + 1 = x+1 → {1, 1}
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr GF49 one_e = GF49::one();
  constexpr auto r = x_e + one_e;
  return r[0] == GF7(1) && r[1] == GF7(1);
}

consteval bool test_gf49_sub()
{
  // (x+1) - x = 1 → {1, 0}
  constexpr GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(1), GF7(1)}}};
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr auto r = xp1_e - x_e;
  return r[0] == GF7(1) && r[1] == GF7(0);
}

consteval bool test_gf49_neg()
{
  // -(x) = {0, 6}  since -(1) ≡ 6 (mod 7)
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr auto r = -x_e;
  return r[0] == GF7(0) && r[1] == GF7(6);
}

consteval bool test_gf49_mul_x_squared()
{
  // x*x = x² ≡ -1 ≡ 6 → {6, 0}
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr auto r = x_e * x_e;
  return r[0] == GF7(6) && r[1] == GF7(0);
}

consteval bool test_gf49_mul_xp1_squared()
{
  // (x+1)² = x²+2x+1 = (-1)+2x+1 = 2x → {0, 2}
  constexpr GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(1), GF7(1)}}};
  constexpr auto r = xp1_e * xp1_e;
  return r[0] == GF7(0) && r[1] == GF7(2);
}

consteval bool test_gf49_mul_zero_absorb()
{
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  return (GF49::zero() * x_e) == GF49::zero();
}

consteval bool test_gf49_mul_one_identity()
{
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  return (GF49::one() * x_e) == x_e;
}

consteval bool test_gf49_inv_x()
{
  // inv(x) = 6x: x*(6x) = 6x² = 6*(-1) ≡ 1 (mod 7) → {0, 6}
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr auto r = inv(x_e);
  return r[0] == GF7(0) && r[1] == GF7(6);
}

consteval bool test_gf49_inv_xp1()
{
  // inv(x+1) = 4+3x: (x+1)(4+3x) = 4+7x+3x² = 4+0+3*(-1) = 1 → {4, 3}
  constexpr GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(1), GF7(1)}}};
  constexpr auto r = inv(xp1_e);
  return r[0] == GF7(4) && r[1] == GF7(3);
}

consteval bool test_gf49_div_self()
{
  // x / x = one()
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  return (x_e / x_e) == GF49::one();
}

consteval bool test_gf49_div_x_by_xp1()
{
  // x / (x+1) = x * inv(x+1) = x*(4+3x) = 4x+3x² = 4x+3*(-1) = 4x-3 ≡ 4x+4 → {4, 4}
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(1), GF7(1)}}};
  constexpr auto r = x_e / xp1_e;
  return r[0] == GF7(4) && r[1] == GF7(4);
}

consteval bool test_gf49_roundtrip()
{
  // (x+2) * inv(x+2) = one()
  constexpr GF49 a{polynomial_nttp<GF7, 1>{{GF7(2), GF7(1)}}};
  return (a * inv(a)) == GF49::one();
}

consteval bool test_gf49_distributivity()
{
  // (x+1) * x should equal 1*x + x*x coefficient-wise
  constexpr GF49 one_e = GF49::one();
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(1), GF7(1)}}};
  return (xp1_e * x_e) == (one_e * x_e + x_e * x_e);
}

consteval bool test_gf49_add_zero_identity()
{
  // x + zero() == x  (additive identity)
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  return (x_e + GF49::zero()) == x_e;
}

consteval bool test_gf49_add_inverse()
{
  // x + (-x) == zero()  (additive inverse)
  constexpr GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  return (x_e + (-x_e)) == GF49::zero();
}

static_assert(test_gf49_zero(), "GF49::zero() must be {0,0}");
static_assert(test_gf49_one(), "GF49::one() must be {1,0}");
static_assert(test_gf49_add_zero_identity(), "GF49: x + zero() == x");
static_assert(test_gf49_add_inverse(), "GF49: x + (-x) == zero()");
static_assert(test_gf49_add(), "GF49: x + 1 == x+1");
static_assert(test_gf49_sub(), "GF49: (x+1) - x == 1");
static_assert(test_gf49_neg(), "GF49: -(x) == 6x");
static_assert(test_gf49_mul_x_squared(), "GF49: x*x == 6");
static_assert(test_gf49_mul_xp1_squared(), "GF49: (x+1)² == 2x");
static_assert(test_gf49_mul_zero_absorb(), "GF49: zero()*x == zero()");
static_assert(test_gf49_mul_one_identity(), "GF49: one()*x == x");
static_assert(test_gf49_inv_x(), "GF49: inv(x) == 6x");
static_assert(test_gf49_inv_xp1(), "GF49: inv(x+1) == 4+3x");
static_assert(test_gf49_div_self(), "GF49: x/x == one()");
static_assert(test_gf49_div_x_by_xp1(), "GF49: x/(x+1) == 4+4x");
static_assert(test_gf49_roundtrip(), "GF49: (x+2)*inv(x+2) == one()");
static_assert(test_gf49_distributivity(), "GF49: distributivity (x+1)*x == x + x²");

// ============================================================
// GF(4) = GF(2)[x]/(x²+x+1)
// Characteristic 2: -x = x, x² ≡ x+1
// ============================================================

consteval bool test_gf4_zero()
{
  constexpr auto z = GF4::zero();
  return z[0] == GF2(0) && z[1] == GF2(0);
}

consteval bool test_gf4_one()
{
  constexpr auto o = GF4::one();
  return o[0] == GF2(1) && o[1] == GF2(0);
}

consteval bool test_gf4_char2_add()
{
  // x + x = 0 in characteristic 2
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  return (x_e + x_e) == GF4::zero();
}

consteval bool test_gf4_char2_neg()
{
  // -x = x in characteristic 2 (since -1 ≡ 1 mod 2)
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  return (-x_e) == x_e;
}

consteval bool test_gf4_mul_x_squared()
{
  // x*x = x² ≡ x+1 mod (x²+x+1) → {1, 1}
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  constexpr GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(1), GF2(1)}}};
  return (x_e * x_e) == xp1_e;
}

consteval bool test_gf4_mul_x_times_xp1()
{
  // x*(x+1) = x²+x ≡ (x+1)+x = 1 in GF(2) → one()
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  constexpr GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(1), GF2(1)}}};
  return (x_e * xp1_e) == GF4::one();
}

consteval bool test_gf4_inv_x()
{
  // inv(x) = x+1: x*(x+1) = 1
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  constexpr GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(1), GF2(1)}}};
  return inv(x_e) == xp1_e;
}

consteval bool test_gf4_inv_xp1()
{
  // inv(x+1) = x: (x+1)*x = 1
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  constexpr GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(1), GF2(1)}}};
  return inv(xp1_e) == x_e;
}

consteval bool test_gf4_div_x_by_xp1()
{
  // x / (x+1) = x * inv(x+1) = x * x = x+1
  constexpr GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(0), GF2(1)}}};
  constexpr GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(1), GF2(1)}}};
  return (x_e / xp1_e) == xp1_e;
}

static_assert(test_gf4_zero(), "GF4::zero() must be {0,0}");
static_assert(test_gf4_one(), "GF4::one() must be {1,0}");
static_assert(test_gf4_char2_add(), "GF4: x + x == zero() (characteristic 2)");
static_assert(test_gf4_char2_neg(), "GF4: -x == x (characteristic 2)");
static_assert(test_gf4_mul_x_squared(), "GF4: x*x == x+1");
static_assert(test_gf4_mul_x_times_xp1(), "GF4: x*(x+1) == one()");
static_assert(test_gf4_inv_x(), "GF4: inv(x) == x+1");
static_assert(test_gf4_inv_xp1(), "GF4: inv(x+1) == x");
static_assert(test_gf4_div_x_by_xp1(), "GF4: x/(x+1) == x+1");

// ============================================================
// GF(289) = GF(17)[x]/(x²+3)
// x² ≡ -3 ≡ 14 (mod 17)
// inv(x) = 11x: 11*14 ≡ 154 ≡ 1 (mod 17)  [154 = 9*17+1]
// inv(x+1) = 13+4x: (x+1)(13+4x) = 13+17x+4x² = 13+0+4*14 = 13+56 ≡ 1 (mod 17) [69=4*17+1]
// inv(3x+5) = 5+14x: (3x+5)(5+14x) = 15x+42x²+25+70x = 25+85x+42*14
//   42*14=588=34*17+10 → 10; 85=5*17 → 0; 25+10=35=2*17+1 → 1 ✓
// ============================================================

consteval bool test_gf289_one()
{
  constexpr auto o = GF289::one();
  return o[0] == GF17(1) && o[1] == GF17(0);
}

consteval bool test_gf289_mul_x_squared()
{
  // x*x = x² ≡ -3 ≡ 14 (mod 17) → {14, 0}
  constexpr GF289 x_e{polynomial_nttp<GF17, 1>{{GF17(0), GF17(1)}}};
  constexpr auto r = x_e * x_e;
  return r[0] == GF17(14) && r[1] == GF17(0);
}

consteval bool test_gf289_inv_x()
{
  // inv(x) = 11x → {0, 11}
  constexpr GF289 x_e{polynomial_nttp<GF17, 1>{{GF17(0), GF17(1)}}};
  constexpr auto r = inv(x_e);
  return r[0] == GF17(0) && r[1] == GF17(11);
}

consteval bool test_gf289_inv_xp1()
{
  // inv(x+1) = 13+4x → {13, 4}
  constexpr GF289 xp1_e{polynomial_nttp<GF17, 1>{{GF17(1), GF17(1)}}};
  constexpr auto r = inv(xp1_e);
  return r[0] == GF17(13) && r[1] == GF17(4);
}

consteval bool test_gf289_inv_3xp5()
{
  // inv(3x+5) = 5+14x → {5, 14}
  constexpr GF289 a{polynomial_nttp<GF17, 1>{{GF17(5), GF17(3)}}};
  constexpr auto r = inv(a);
  return r[0] == GF17(5) && r[1] == GF17(14);
}

consteval bool test_gf289_roundtrip()
{
  // (x+1) * inv(x+1) = one()
  constexpr GF289 xp1_e{polynomial_nttp<GF17, 1>{{GF17(1), GF17(1)}}};
  return (xp1_e * inv(xp1_e)) == GF289::one();
}

static_assert(test_gf289_one(), "GF289::one() must be {1,0}");
static_assert(test_gf289_mul_x_squared(), "GF289: x*x == 14");
static_assert(test_gf289_inv_x(), "GF289: inv(x) == 11x");
static_assert(test_gf289_inv_xp1(), "GF289: inv(x+1) == 13+4x");
static_assert(test_gf289_inv_3xp5(), "GF289: inv(3x+5) == 5+14x");
static_assert(test_gf289_roundtrip(), "GF289: (x+1)*inv(x+1) == one()");

// ============================================================
// GF(25) = GF(5)[x]/(x²+2)
// x² ≡ -2 ≡ 3 (mod 5)
// inv(x) = 2x: x*(2x) = 2x² = 2*3 = 6 ≡ 1 (mod 5) ✓
// ============================================================

consteval bool test_gf25_one()
{
  constexpr auto o = GF25::one();
  return o[0] == GF5(1) && o[1] == GF5(0);
}

consteval bool test_gf25_mul_x_squared()
{
  // x*x = x² ≡ -2 ≡ 3 (mod 5) → {3, 0}
  constexpr GF25 x_e{polynomial_nttp<GF5, 1>{{GF5(0), GF5(1)}}};
  constexpr auto r = x_e * x_e;
  return r[0] == GF5(3) && r[1] == GF5(0);
}

consteval bool test_gf25_inv_x()
{
  // inv(x) = 2x → {0, 2}
  constexpr GF25 x_e{polynomial_nttp<GF5, 1>{{GF5(0), GF5(1)}}};
  constexpr auto r = inv(x_e);
  return r[0] == GF5(0) && r[1] == GF5(2);
}

consteval bool test_gf25_roundtrip()
{
  // x / x = one()
  constexpr GF25 x_e{polynomial_nttp<GF5, 1>{{GF5(0), GF5(1)}}};
  return (x_e / x_e) == GF25::one();
}

static_assert(test_gf25_one(), "GF25::one() must be {1,0}");
static_assert(test_gf25_mul_x_squared(), "GF25: x*x == 3");
static_assert(test_gf25_inv_x(), "GF25: inv(x) == 2x");
static_assert(test_gf25_roundtrip(), "GF25: x/x == one()");

// ============================================================
// Runtime tests (volatile-sourced input to defeat compile-time eval)
// ============================================================

bool runtime_test_gf49_add()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(v0), GF7(v1)}}};
  GF49 one_e = GF49::one();
  auto r = x_e + one_e;
  return r[0] == GF7(1) && r[1] == GF7(1);
}

bool runtime_test_gf49_mul_x_squared()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(v0), GF7(v1)}}};
  auto r = x_e * x_e;
  return r[0] == GF7(6) && r[1] == GF7(0);
}

bool runtime_test_gf49_inv_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(v0), GF7(v1)}}};
  auto r = inv(x_e);
  return r[0] == GF7(0) && r[1] == GF7(6);
}

bool runtime_test_gf49_div()
{
  volatile std::uint64_t v1 = 1;
  GF49 x_e{polynomial_nttp<GF7, 1>{{GF7(0), GF7(v1)}}};
  GF49 xp1_e{polynomial_nttp<GF7, 1>{{GF7(v1), GF7(v1)}}};
  auto r = x_e / xp1_e;
  return r[0] == GF7(4) && r[1] == GF7(4);
}

bool runtime_test_gf4_inv_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  GF4 x_e{polynomial_nttp<GF2, 1>{{GF2(v0), GF2(v1)}}};
  GF4 xp1_e{polynomial_nttp<GF2, 1>{{GF2(v1), GF2(v1)}}};
  return inv(x_e) == xp1_e;
}

bool runtime_test_gf289_inv_x()
{
  volatile std::uint64_t v0 = 0;
  volatile std::uint64_t v1 = 1;
  GF289 x_e{polynomial_nttp<GF17, 1>{{GF17(v0), GF17(v1)}}};
  auto r = inv(x_e);
  return r[0] == GF17(0) && r[1] == GF17(11);
}

// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

int main()
{
  // Collect compile-time results
  constexpr bool ct_gf49_zero = test_gf49_zero();
  constexpr bool ct_gf49_one = test_gf49_one();
  constexpr bool ct_gf49_add_zero_id = test_gf49_add_zero_identity();
  constexpr bool ct_gf49_add_inv = test_gf49_add_inverse();
  constexpr bool ct_gf49_add = test_gf49_add();
  constexpr bool ct_gf49_sub = test_gf49_sub();
  constexpr bool ct_gf49_neg = test_gf49_neg();
  constexpr bool ct_gf49_mul_xx = test_gf49_mul_x_squared();
  constexpr bool ct_gf49_mul_xp1_sq = test_gf49_mul_xp1_squared();
  constexpr bool ct_gf49_zero_absorb = test_gf49_mul_zero_absorb();
  constexpr bool ct_gf49_one_id = test_gf49_mul_one_identity();
  constexpr bool ct_gf49_inv_x = test_gf49_inv_x();
  constexpr bool ct_gf49_inv_xp1 = test_gf49_inv_xp1();
  constexpr bool ct_gf49_div_self = test_gf49_div_self();
  constexpr bool ct_gf49_div = test_gf49_div_x_by_xp1();
  constexpr bool ct_gf49_roundtrip = test_gf49_roundtrip();
  constexpr bool ct_gf49_distrib = test_gf49_distributivity();

  constexpr bool ct_gf4_zero = test_gf4_zero();
  constexpr bool ct_gf4_one = test_gf4_one();
  constexpr bool ct_gf4_char2_add = test_gf4_char2_add();
  constexpr bool ct_gf4_char2_neg = test_gf4_char2_neg();
  constexpr bool ct_gf4_mul_xx = test_gf4_mul_x_squared();
  constexpr bool ct_gf4_mul_x_xp1 = test_gf4_mul_x_times_xp1();
  constexpr bool ct_gf4_inv_x = test_gf4_inv_x();
  constexpr bool ct_gf4_inv_xp1 = test_gf4_inv_xp1();
  constexpr bool ct_gf4_div = test_gf4_div_x_by_xp1();

  constexpr bool ct_gf289_one = test_gf289_one();
  constexpr bool ct_gf289_mul_xx = test_gf289_mul_x_squared();
  constexpr bool ct_gf289_inv_x = test_gf289_inv_x();
  constexpr bool ct_gf289_inv_xp1 = test_gf289_inv_xp1();
  constexpr bool ct_gf289_inv_3xp5 = test_gf289_inv_3xp5();
  constexpr bool ct_gf289_roundtrip = test_gf289_roundtrip();

  constexpr bool ct_gf25_one = test_gf25_one();
  constexpr bool ct_gf25_mul_xx = test_gf25_mul_x_squared();
  constexpr bool ct_gf25_inv_x = test_gf25_inv_x();
  constexpr bool ct_gf25_roundtrip = test_gf25_roundtrip();

  // Runtime results
  bool rt_gf49_add = runtime_test_gf49_add();
  bool rt_gf49_mul_xx = runtime_test_gf49_mul_x_squared();
  bool rt_gf49_inv_x = runtime_test_gf49_inv_x();
  bool rt_gf49_div = runtime_test_gf49_div();
  bool rt_gf4_inv_x = runtime_test_gf4_inv_x();
  bool rt_gf289_inv_x = runtime_test_gf289_inv_x();

  constexpr bool all_ct =
    ct_gf49_zero && ct_gf49_one && ct_gf49_add_zero_id && ct_gf49_add_inv && ct_gf49_add && ct_gf49_sub &&
    ct_gf49_neg && ct_gf49_mul_xx && ct_gf49_mul_xp1_sq && ct_gf49_zero_absorb && ct_gf49_one_id && ct_gf49_inv_x &&
    ct_gf49_inv_xp1 && ct_gf49_div_self && ct_gf49_div && ct_gf49_roundtrip && ct_gf49_distrib && ct_gf4_zero &&
    ct_gf4_one && ct_gf4_char2_add && ct_gf4_char2_neg && ct_gf4_mul_xx && ct_gf4_mul_x_xp1 && ct_gf4_inv_x &&
    ct_gf4_inv_xp1 && ct_gf4_div && ct_gf289_one && ct_gf289_mul_xx && ct_gf289_inv_x && ct_gf289_inv_xp1 &&
    ct_gf289_inv_3xp5 && ct_gf289_roundtrip && ct_gf25_one && ct_gf25_mul_xx && ct_gf25_inv_x && ct_gf25_roundtrip;

  if constexpr (all_ct)
    std::println("All compile-time tests passed.");

  bool all_rt = rt_gf49_add && rt_gf49_mul_xx && rt_gf49_inv_x && rt_gf49_div && rt_gf4_inv_x && rt_gf289_inv_x;

  if (all_rt)
    std::println("All runtime tests passed.");

  return (all_ct && all_rt) ? 0 : 1;
}
