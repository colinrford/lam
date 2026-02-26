/*
 *  finite_field_traits_test.cpp – written by Colin Ford
 *    see github.com/colinrford/lam for license info
 *
 *  Tests for finite_field_traits specializations in lam-interop:
 *    - ZqElement<T, P>                           (regression: field_order rename)
 *    - finite_field_extension<K, N, Irreducible> (Phase 3)
 *
 *  Verified properties:
 *    is_finite_field == true
 *    field_order     == expected value (p for GF(p), p^N for GF(p^N))
 *    mul/add/sub     delegate correctly to operators
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;

using namespace lam::cbn;

using GF2 = ZqElement<std::uint64_t, 2>;
using GF5 = ZqElement<std::uint64_t, 5>;
using GF7 = ZqElement<std::uint64_t, 7>;
using GF17 = ZqElement<std::uint64_t, 17>;

using lam::polynomial::univariate::finite_field_traits;

// ============================================================
// ZqElement regression: field_order rename must not break existing specialization
// ============================================================

static_assert(finite_field_traits<GF7>::is_finite_field, "GF7 must be a finite field");
static_assert(finite_field_traits<GF7>::field_order == 7, "GF7 field_order must be 7");
static_assert(finite_field_traits<GF17>::is_finite_field, "GF17 must be a finite field");
static_assert(finite_field_traits<GF17>::field_order == 17, "GF17 field_order must be 17");

// ============================================================
// finite_field_extension specialization
// ============================================================

constexpr lam::polynomial_nttp<GF2, 2> gf4_mod{{GF2(1), GF2(1), GF2(1)}};
constexpr lam::polynomial_nttp<GF7, 2> gf49_mod{{GF7(1), GF7(0), GF7(1)}};
constexpr lam::polynomial_nttp<GF5, 2> gf25_mod{{GF5(2), GF5(0), GF5(1)}};
constexpr lam::polynomial_nttp<GF17, 2> gf289_mod{{GF17(3), GF17(0), GF17(1)}};
constexpr lam::polynomial_nttp<GF7, 3> gf343_mod{{GF7(1), GF7(1), GF7(0), GF7(1)}}; // x³+x+1

using GF4 = lam::finite_field_extension<GF2, 2, gf4_mod>;
using GF49 = lam::finite_field_extension<GF7, 2, gf49_mod>;
using GF25 = lam::finite_field_extension<GF5, 2, gf25_mod>;
using GF289 = lam::finite_field_extension<GF17, 2, gf289_mod>;
using GF343 = lam::finite_field_extension<GF7, 3, gf343_mod>; // degree-3, exercises loop

// is_finite_field
static_assert(finite_field_traits<GF4>::is_finite_field, "GF4 must be a finite field");
static_assert(finite_field_traits<GF49>::is_finite_field, "GF49 must be a finite field");
static_assert(finite_field_traits<GF25>::is_finite_field, "GF25 must be a finite field");
static_assert(finite_field_traits<GF289>::is_finite_field, "GF289 must be a finite field");
static_assert(finite_field_traits<GF343>::is_finite_field, "GF343 must be a finite field");

// field_order = p^N — degree-2 cases and one degree-3 to exercise the loop fully
static_assert(finite_field_traits<GF4>::field_order == 4, "GF4 field_order must be 2^2 = 4");
static_assert(finite_field_traits<GF49>::field_order == 49, "GF49 field_order must be 7^2 = 49");
static_assert(finite_field_traits<GF25>::field_order == 25, "GF25 field_order must be 5^2 = 25");
static_assert(finite_field_traits<GF289>::field_order == 289, "GF289 field_order must be 17^2 = 289");
static_assert(finite_field_traits<GF343>::field_order == 343, "GF343 field_order must be 7^3 = 343");

// TODO: test field_order = 0 fallback when K has no finite_field_traits specialization
// (requires a dummy field_element_c_weak type with no specialization; deferred)

// mul/add/sub delegate to operators — spot-check one ring
consteval bool test_traits_arithmetic_gf49()
{
  constexpr GF49 x_e{lam::polynomial_nttp<GF7, 1>{{GF7(0), GF7(1)}}};
  constexpr GF49 one_e = GF49::one();
  using T = finite_field_traits<GF49>;
  constexpr auto sum = T::add(x_e, one_e);  // x + 1 = {1, 1}
  constexpr auto diff = T::sub(x_e, one_e); // x - 1 = {6, 1} in GF(7)
  constexpr auto prod = T::mul(x_e, x_e);   // x*x = {6, 0}  (x² ≡ -1 ≡ 6)
  return sum[0] == GF7(1) && sum[1] == GF7(1) && diff[0] == GF7(6) && diff[1] == GF7(1) && prod[0] == GF7(6) &&
         prod[1] == GF7(0);
}
static_assert(test_traits_arithmetic_gf49(), "finite_field_traits<GF49> mul/add/sub must match operators");

int main()
{
  std::println("All finite_field_traits static assertions passed.");
  return 0;
}
