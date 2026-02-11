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
  static constexpr std::size_t modulus = static_cast<std::size_t>(P);

  static constexpr K mul(const K& a, const K& b) {
      if constexpr (P == 0xFFFFFFFF00000001ULL) {
          // Fast Solinas path
          unsigned __int128 prod = static_cast<unsigned __int128>(a.data[0]) * b.data[0];
          return K(lam::polynomial::univariate::ntt::reduce_solinas(prod));
      } else {
          return a * b;
      }
  }
};
