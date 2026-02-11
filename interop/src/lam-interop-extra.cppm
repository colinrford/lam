/*
 *  lam-interop-extra.cppm
 *    Advanced interoperability utilities
 */

export module lam.interop.extra;

import std;
import lam.concepts;
import lam.polynomial_nttp;
import lam.ctbignum;

namespace lam::interop
{

/**
 *  Helper to create a polynomial from a vector of string representations of coefficients.
 *  Useful for constructing large finite field polynomials from text data.
 */
export template<typename Zq, std::size_t N>
constexpr auto make_poly_from_strings(const std::vector<std::string>& coeffs) -> lam::polynomial::univariate::polynomial_nttp<Zq, N>
{
    lam::polynomial::univariate::polynomial_nttp<Zq, N> p{};
    for (std::size_t i = 0; i < coeffs.size() && i <= N; ++i)
    {
        // Parse string to uint64_t or whatever the underlying type is (assuming ZqElement)
        // This is a simple implementation for common ZqElement<uint64_t>
        // For full generic big int support, we'd need more complex parsing logic from ctbignum if exposed.
        // Assuming ctbignum elements can be constructed from integer literals.
        
        // Using stoull for basic 64-bit primes.
        // If ZqElement supports string parsing directly, we should use that.
        // ctbignum literals are _Z, so let's rely on constructing from raw integers number.
        
        std::uint64_t val = std::stoull(coeffs[i]);
        p.coefficients[i] = Zq(val); 
    }
    return p;
}

/**
 *  Print a polynomial with finite field coefficients in a readable format.
 */
export template<typename Zq, std::size_t N>
void print_poly(const lam::polynomial::univariate::polynomial_nttp<Zq, N>& p, std::string_view name)
{
    std::print("{} = ", name);
    bool first = true;
    for (std::size_t i = 0; i <= N; ++i)
    {
        // Check if zero (assuming .data internal access or relying on operator==)
        // ZqElement usually has .data for the raw value
        auto val = p[i];
        // We can cast to its integer type if we know it? 
        // Or just rely on format if available.
        // ctbignum elements usually format as hex.
        
        // Let's assume we want to print non-zero terms
        // But for finite fields, "zero" comparison might be tricky without a clear "is_zero".
        // Let's just print everything for small N or iterate.
        
        // To be safe and generic, let's strict to simple printing.
         if (i == 0) std::print("{}", (std::uint64_t)val.data[0]); // Hacky access to data[0] assuming typical ctbignum limb 
         else std::print(" + {}*x^{}", (std::uint64_t)val.data[0], i);
    }
    std::println("");
}

}
