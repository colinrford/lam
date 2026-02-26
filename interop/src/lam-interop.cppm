/*
 *  lam-interop.cppm
 *    see github.com/colinrford/lam for license info
 *
 *  Interoperability module for lam libraries
 */

export module lam.interop;

import std;
import lam.concepts;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.bitvector;

template<typename T, auto P>
struct lam::polynomial::univariate::finite_field_traits<lam::cbn::ZqElement<T, P>>
{
  using K = lam::cbn::ZqElement<T, P>;
  static constexpr bool is_finite_field = true;
  static constexpr std::size_t field_order = static_cast<std::size_t>(P);

  static constexpr K mul(const K& a, const K& b) {
      if constexpr (P == 0xFFFFFFFF00000001ULL) {
          // Fast Solinas path
          unsigned __int128 prod = static_cast<unsigned __int128>(a.data[0]) * b.data[0];
          return K(lam::polynomial::univariate::ntt::reduce_solinas(prod));
      } else {
          return a * b;
      }
  }
  static constexpr K add(const K& a, const K& b) { return a + b; }
  static constexpr K sub(const K& a, const K& b) { return a - b; }
};

// finite_field_traits specialization for finite_field_extension<K, N, Irreducible>.
// field_order is |K|^N, derived from finite_field_traits<K>::field_order when K is a finite field.
// when K itself is a finite field (e.g. ZqElement). Defaults to 0 otherwise.
template<lam::concepts::experimental::field_element_c_weak K,
         std::size_t N,
         lam::polynomial_nttp<K, N> Irreducible>
struct lam::polynomial::univariate::finite_field_traits<lam::finite_field_extension<K, N, Irreducible>>
{
  using F = lam::finite_field_extension<K, N, Irreducible>;

  static constexpr bool is_finite_field = true;

  static constexpr std::size_t field_order = []() -> std::size_t {
    using KTraits = lam::polynomial::univariate::finite_field_traits<K>;
    if constexpr (KTraits::is_finite_field) {
      std::size_t result = 1;
      for (std::size_t i = 0; i < N; ++i)
        result *= KTraits::field_order;
      return result;
    }
    return 0;
  }();

  static constexpr F mul(const F& a, const F& b) { return a * b; }
  static constexpr F add(const F& a, const F& b) { return a + b; }
  static constexpr F sub(const F& a, const F& b) { return a - b; }
};
