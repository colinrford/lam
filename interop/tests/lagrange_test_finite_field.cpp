/*
 *  lagrange_test_finite_field.cpp
 *  unit test for lagrange interpolation over finite fields
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam;
using namespace lam::cbn;
using namespace lam::cbn::literals;

int main()
{
  int failures = 0;

  // Define a small field GF(17)
  using GF17 = ZqElement<std::uint64_t, 17>;

  // ============================================================
  // Test 1: Linear in GF(17) -> y = x
  // Points: (0,0), (1,1)
  // ============================================================
  {
    std::array<GF17, 2> x = {GF17(0), GF17(1)};
    std::array<GF17, 2> y = {GF17(0), GF17(1)};
    auto poly = lagrange_interpolate<GF17, 1>(x, y);

    // Expected: 0 + 1*x
    if (poly.coefficients[0] != GF17(0) || poly.coefficients[1] != GF17(1))
    {
      std::println("FAIL: GF(17) linear");
      failures++;
    }
  }

  // ============================================================
  // Test 2: Quadratic in GF(17) -> y = x^2 + 1
  // Points: (0, 1), (1, 2), (2, 5)
  // 0 -> 0^2 + 1 = 1
  // 1 -> 1^2 + 1 = 2
  // 2 -> 2^2 + 1 = 5
  // ============================================================
  {
    std::array<GF17, 3> x = {GF17(0), GF17(1), GF17(2)};
    std::array<GF17, 3> y = {GF17(1), GF17(2), GF17(5)};
    auto poly = lagrange_interpolate<GF17, 2>(x, y);

    // Expected: 1 + 0*x + 1*x^2 -> [1, 0, 1]
    if (poly.coefficients[0] != GF17(1) || poly.coefficients[1] != GF17(0) || poly.coefficients[2] != GF17(1))
    {
      std::println("FAIL: GF(17) quadratic");
      failures++;
    }
  }

  // ============================================================
  // Test 3: Collinear points for Quadratic (Effectively Linear) in GF(17)
  // Points: (0, 1), (1, 2), (2, 3)
  // ============================================================
  {
    std::array<GF17, 3> x = {GF17(0), GF17(1), GF17(2)};
    std::array<GF17, 3> y = {GF17(1), GF17(2), GF17(3)};
    auto poly = lagrange_interpolate<GF17, 2>(x, y);

    // Expected: 1 + 1*x + 0*x^2 -> [1, 1, 0]
    if (poly.coefficients[0] != GF17(1)
     || poly.coefficients[1] != GF17(1)
     || poly.coefficients[2] != GF17(0))
    {
      std::println("FAIL: GF(17) effectively linear");
      failures++;
    }
  }

  // ============================================================
  // Test 4: Compile-time evaluation in GF(17)
  // Points: (0, 0), (1, 1) -> y = x
  // ============================================================
  {
    constexpr std::array<GF17, 2> x_c = {GF17(0), GF17(1)};
    constexpr std::array<GF17, 2> y_c = {GF17(0), GF17(1)};
    constexpr auto poly_c = lagrange_interpolate<GF17, 1>(x_c, y_c);

    static_assert(poly_c.degree == 1, "Degree must be 1");
    static_assert(poly_c.coefficients[0] == GF17(0), "Coefficient 0 must be 0");
    static_assert(poly_c.coefficients[1] == GF17(1), "Coefficient 1 must be 1");
  }

  // ============================================================
  // Test 5: Exactness Stress Test in GF(4179340454199820289)
  // Large prime, Degree 10 polynomial recovery
  // ============================================================
  {
    constexpr auto large_prime = 4179340454199820289_Z;
    using GFLarge = decltype(lam::cbn::Zq(large_prime));

    // P(x) = 1 + 2x + 3x^2 + ... + 11x^10
    auto expected_poly = [&]() {
      polynomial_nttp<GFLarge, 10> p;
      for (std::size_t i = 0; i <= 10; ++i)
      {
        p.coefficients[i] = GFLarge(i + 1);
      }
      return p;
    }();

    std::array<GFLarge, 11> x_vals;
    std::array<GFLarge, 11> y_vals;
    for (std::size_t i = 0; i < 11; ++i)
    {
      // 0, 1, 4, 9, 16 ...
      x_vals[i] = GFLarge(i * i);
    }
    // ensure unique points
    for (std::size_t i = 0; i < 11; ++i)
      x_vals[i] = GFLarge(i * 3 + 1);

    for (std::size_t i = 0; i < 11; ++i)
    {
      y_vals[i] = expected_poly(x_vals[i]);
    }

    auto recovered = lagrange_interpolate<GFLarge, 10>(x_vals, y_vals);

    bool exact = true;
    for (std::size_t i = 0; i <= 10; ++i)
    {
      if (recovered.coefficients[i] != expected_poly.coefficients[i])
      {
        exact = false;
        break;
      }
    }

    if (!exact)
    {
      std::println("FAIL: GF(Large) degree-10 recovery");
      failures++;
    }
  }

  // ============================================================
  // Test 6: Max Degree Interpolation in GF(17)
  // Degree 16 interpolation utilizing all 17 points
  // ============================================================
  {
    std::array<GF17, 17> x_vals;
    std::array<GF17, 17> y_vals;
    for (std::size_t i = 0; i < 17; ++i)
    {
      x_vals[i] = GF17(i);
      // Let's interpolate y = x^16. By Fermat's Little Theorem, x^16 = 1 for x
      // != 0, and 0 for x = 0.
      y_vals[i] = i == 0 ? GF17(0) : GF17(1);
    }

    auto recovered = lagrange_interpolate<GF17, 16>(x_vals, y_vals);

    // In GF(17), x^16 = 1 (x!=0), 0 (x=0).
    // The polynomial that is 0 at 0 and 1 elsewhere is exactly P(x) = x^16.
    // So coefficients should be 0 except for x^16 which is 1.
    bool exact = true;
    for (std::size_t i = 0; i < 16; ++i)
    {
      if (recovered.coefficients[i] != GF17(0))
        exact = false;
    }
    if (recovered.coefficients[16] != GF17(1))
      exact = false;

    if (!exact)
    {
      std::println("FAIL: GF(17) max-degree (16) recovery");
      failures++;
    }
  }

  return failures;
}
