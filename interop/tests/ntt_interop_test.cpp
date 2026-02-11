/*
 *  ntt_interop_test.cpp
 *  Verifies that importing lam.interop correctly provides traits for ctbignum
 * compatibility.
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;

using namespace lam::cbn::literals;

// Prime: 29 * 2^57 + 1
constexpr auto ntt_prime = 4179340454199820289_Z;
using Field = decltype(lam::cbn::Zq(ntt_prime));

using namespace lam::polynomial::univariate::ntt;

int main() {
  std::println("Testing lam.interop integration...");

  std::vector<Field> data = {Field(1), Field(2), Field(3), Field(4)};
  auto original = data;

  // Should succeed if interop module works
  ntt_transform(data, false);
  ntt_transform(data, true);

  for (std::size_t i = 0; i < 4; ++i) {
    if (data[i] != original[i]) {
      std::println("Interop Round-Trip Failed!");
      return 1;
    }
  }

  std::println("Interop Integration Verified.");

  // 2. High-Level Polynomial Multiplication Test
  std::println("Testing High-Level Polynomial Multiplication (a * b)...");

  // P = 4179340454199820289
  // a = x + 1
  // b = x - 1
  // a * b = x^2 - 1 = x^2 + (P - 1)

  lam::polynomial::polynomial_nttp<Field, 64> a;
  lam::polynomial::polynomial_nttp<Field, 64> b;

  a.coefficients[1] = Field(1);
  a.coefficients[0] = Field(1); // x + 1
  // x - 1 equivalent to x + (P - 1)
  b.coefficients[1] = Field(1);
  b.coefficients[0] = Field(4179340454199820288ULL);

  auto product = a * b;

  // Expected: coefficient[2] = 1, coefficient[0] = -1 (P-1)
  if (product.coefficients[2] != Field(1)) {
    std::println("Polynomial Mul Failed: x^2 coeff wrong");
    return 1;
  }
  if (product.coefficients[0] != Field(4179340454199820288ULL)) {
    std::println("Polynomial Mul Failed: constant coeff wrong");
    return 1;
  }

  std::println("Polynomial Multiplication Verified (NTT Triggered).");

  // 3. Randomized Stress Test vs Naive Multiplication
  std::println("Running Randomized Stress Test (N=128)...");

  constexpr std::size_t N = 128;
  lam::polynomial::polynomial_nttp<Field, N> p1;
  lam::polynomial::polynomial_nttp<Field, N> p2;

  // Deterministic pseudo-random generation for reproducibility
  auto lcg = [state = 12345ULL]() mutable {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return state;
  };

  for (std::size_t i = 0; i <= N; ++i)
    p1.coefficients[i] = Field(lcg());
  for (std::size_t i = 0; i <= N; ++i)
    p2.coefficients[i] = Field(lcg());

  // Fast Multiplication (NTT)
  auto p_ntt = p1 * p2;

  // Naive Multiplication (Reference)
  lam::polynomial::polynomial_nttp<Field, 2 * N> p_naive{};
  for (std::size_t i = 0; i <= N; ++i) {
    for (std::size_t j = 0; j <= N; ++j) {
      p_naive.coefficients[i + j] = p_naive.coefficients[i + j] +
                                    (p1.coefficients[i] * p2.coefficients[j]);
    }
  }

  // Verify
  for (std::size_t i = 0; i <= 2 * N; ++i) {
    if (p_ntt.coefficients[i] != p_naive.coefficients[i]) {
      std::println("Stress Test Failed at index {}!", i);
      // std::println("Expected: {}", p_naive.coefficients[i]);
      // std::println("Got:      {}", p_ntt.coefficients[i]);
      return 1;
    }
  }

  std::println("Randomized Stress Test Passed (N=128).");
  return 0;
}
