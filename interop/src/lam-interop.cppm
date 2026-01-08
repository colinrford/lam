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
  static constexpr bool is_finite_field = true;
  static constexpr std::size_t modulus = static_cast<std::size_t>(P);
};
