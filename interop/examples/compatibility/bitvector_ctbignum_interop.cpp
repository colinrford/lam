/*
 *  ctbignum_interop.cpp – written by Colin Ford
 *    see github.com/colinrford/bitvector
 *    currently unlicensed
 */

import std;
import lam.bitvector;
import lam.ctbignum;

using namespace lam::bitvec;
using namespace lam::cbn;

void test_big_int_roundtrip() {
  std::println("Testing big_int <-> bitvector roundtrip...");

  constexpr std::uint64_t TEST_VAL_1 = 0xDEADBEEFULL;
  constexpr std::uint64_t TEST_VAL_2 = 0xCAFEBABEULL;
  constexpr std::uint64_t TEST_VAL_3 = 0x123456789ABCDEF0ULL;
  constexpr std::uint64_t TEST_VAL_4 = 0xFEDCBA9876543210ULL;
  big_int<4> original{TEST_VAL_1, TEST_VAL_2, TEST_VAL_3, TEST_VAL_4};
  bitvector bv(original);

  constexpr std::size_t EXPECTED_BITS = 4ULL * 64;
  if (bv.size() != EXPECTED_BITS)
    throw std::runtime_error(std::format("Size mismatch: expected {}, got {}",
                                         EXPECTED_BITS, bv.size()));

  big_int<4> recovered{};
  bv.export_words(recovered);

  if (original != recovered) {
    std::println("Mismatch detected!");
    std::println("Original:  {:016x} {:016x} {:016x} {:016x}", original[0],
                 original[1], original[2], original[3]);
    std::println("Recovered: {:016x} {:016x} {:016x} {:016x}", recovered[0],
                 recovered[1], recovered[2], recovered[3]);
    throw std::runtime_error("Data mismatch in big_int roundtrip");
  }

  std::println("big_int roundtrip passed.");
}

void test_zq_element_conversion() {
  std::println("Testing ZqElement -> bitvector conversion...");

  constexpr std::uint64_t ZQ_MODULUS = 0xFFFFFFFFFFFFFFC5ULL;
  using MyZq = ZqElement<std::uint64_t, ZQ_MODULUS>;
  constexpr std::uint64_t ZQ_TEST_VAL = 123456789ULL;
  MyZq zq_val(ZQ_TEST_VAL);

  // ZqElement usually has an explicit operator to get the big_int data or just
  // exposes .data Based on user snippet: explicit operator auto() const {
  // return data; } We can convert to big_int first, OR if ZqElement exposes a
  // range interface (it has .data member which is big_int) user snippet said:
  // "big_int<sizeof...(Modulus), T> data;"

  bitvector bv(zq_val.data);

  constexpr std::size_t EXPECTED_BITS_ZQ = 1ULL * 64;
  if (bv.size() != EXPECTED_BITS_ZQ)
    throw std::runtime_error("Size mismatch for ZqElement");

  big_int<1> recovered_data{};
  bv.export_words(recovered_data);

  if (recovered_data != zq_val.data)
    throw std::runtime_error("Data mismatch in ZqElement conversion");

  std::println("ZqElement conversion passed.");
}

consteval bool test_constexpr_interop() {
  constexpr std::uint64_t CONSTEXPR_TEST_VAL_1 = 0xAAAA5555AAAA5555ULL;
  constexpr std::uint64_t CONSTEXPR_TEST_VAL_2 = 0xBBBB6666BBBB6666ULL;
  big_int<2> original{CONSTEXPR_TEST_VAL_1, CONSTEXPR_TEST_VAL_2};
  bitvector bv(original);

  constexpr std::size_t EXPECTED_CONSTEXPR_BITS = 128;
  if (bv.size() != EXPECTED_CONSTEXPR_BITS)
    return false;

  const auto *words = bv.data();
  if (words[0] != original[0])
    return false;
  if (words[1] != original[1])
    return false;

  for (std::size_t i = 0; i < EXPECTED_CONSTEXPR_BITS; ++i) {
    constexpr std::size_t BITS_PER_WORD = 64;
    std::size_t word_idx = i / BITS_PER_WORD;
    std::size_t bit_idx = i % BITS_PER_WORD;

    bool expected = (original[word_idx] >> bit_idx) & 1;
    if (bv[i] != expected)
      return false;
  }

  big_int<2> recovered{};
  bv.export_words(recovered);

  for (int i = 0; i < 2; ++i)
    if (original[i] != recovered[i])
      return false;

  return true;
}

// Verify compile-time execution
static_assert(test_constexpr_interop(), "Constexpr interop failed!");

int main() {
  try {
    test_big_int_roundtrip();
    test_zq_element_conversion();
  } catch (const std::exception &e) {
    std::println("Test failed: {}", e.what());
    return 1;
  }

  std::println("All interoperability tests passed.");
  return 0;
}
