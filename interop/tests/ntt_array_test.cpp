/*
 *  ntt_array_test.cpp
 *  Verifies that NTT works with std::array (requires generic range support).
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;

using namespace lam::cbn::literals;

constexpr auto ntt_prime = 4179340454199820289_Z;
using Field = decltype(lam::cbn::Zq(ntt_prime));

using namespace lam::polynomial::univariate::ntt;

int main()
{
  // Test with std::array (N=4)
  std::array<Field, 4> data = {Field(1), Field(2), Field(3), Field(4)};
  auto original = data;

  // This check should pass if the interface accepts generic ranges
  ntt_transform(data, false);
  ntt_transform(data, true);

  for (std::size_t i = 0; i < 4; ++i)
  {
    if (data[i] != original[i])
    {
      std::println("Array round trip failed!");
      return 1;
    }
  }
  std::println("Array NTT passed.");
  return 0;
}
