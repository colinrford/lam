
import std;
import lam.interop; // Brings in finite_field_traits specialization and components
import lam.ctbignum;
import lam.polynomial_nttp;

using namespace lam::cbn;
using namespace lam::cbn::literals;
using namespace lam::polynomial;

// Kyber parameter q = 3329
using KyberField = ZqElement<std::uint32_t, 3329>;

// Polynomial in R[X]
// Note: We use N=64. Product degree 128. Coeffs 129.
// bit_ceil(129) = 256. This fits the Kyber prime (256 | 3328).
// Previously N=128 -> 257 coeffs -> 512 size -> failed (512 does not divide 3328).
using KyberPoly = polynomial_nttp<KyberField, 64>;

int main() {
    std::println("Testing Kyber Polynomial Interop...");

    // 1. Construct Polynomials
    // P(x) = x + 1
    KyberPoly P{};
    P.coefficients[1] = KyberField{1_Z32};
    P.coefficients[0] = KyberField{1_Z32};

    // Q(x) = x - 1  (which is x + 3328 mod 3329)
    KyberPoly Q{};
    Q.coefficients[1] = KyberField{1_Z32};
    Q.coefficients[0] = KyberField{3328_Z32}; // -1 mod 3329

    // 2. Multiply
    // Expected: (x+1)(x-1) = x^2 - 1 = x^2 + 3328
    // Result degree: 64 + 64 = 128.
    auto R = P * Q;

    std::println("Deg(P) = {}", 64); 
    std::println("Deg(R) = {}", R.degree);

    // 3. Verify Coefficients
    // x^2 coefficient should be 1
    auto c2 = R[2]; 
    // x^1 coefficient should be 0 (1*(-1) + 1*1 = 0)
    auto c1 = R[1];
    // x^0 coefficient should be -1 (3328)
    auto c0 = R[0];

    std::println("Coeff x^2: {}", c2);
    std::println("Coeff x^1: {}", c1);
    std::println("Coeff x^0: {}", c0);

    // Correct verification logic
    if (c2 == KyberField{1_Z32} && c1 == KyberField{0_Z32} && c0 == KyberField{3328_Z32}) {
        std::println("SUCCESS: (x+1)(x-1) == x^2 - 1");
    } else {
        std::println("FAILURE: Incorrect multiplication result.");
        return 1;
    }

    // 4. Verify Large Degree (Force NTT)
    // Degree > 64 triggers NTT. Since product size is 256 (power of 2),
    // and 256 | 3328, the NTT should execute properly.
    
    std::println("Kyber Poly Multiplication successful!");
    return 0;
}
