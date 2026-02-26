# GF extension fields — plan (lam work)

Option B needs three things in lam, in order:

- 1. `poly_rem<Mod>(a)` in `polynomial_nttp`
  Synthetic division with a compile-time modulus. Numerator is runtime 
  (`polynomial_nttp<K, M>`), denominator is NTTP (Irreducible). Returns 
  remainder of degree < N. Pure coefficient arithmetic —
  no FFT/NTT involved, just a loop. Leading coefficient of Mod is known at 
  compile time so the division is exact.

- 2. `gf_ext_element<K, Irreducible>` struct
  Wraps `polynomial_nttp<K, N-1>` (the representative, degree < N). 
  Operators:
  - +, -: coefficient-wise, already works through `polynomial_nttp`
  - \*: polynomial multiply → `poly_rem<Irreducible>` to reduce
  - inv: extended Euclidean in K[x] mod Irreducible

  The irreducible being a NTTP means extension field structure is fully 
  compile-time. Same trade-off as `polynomial_nttp` degree — need a dispatch 
  bridge in lambot.

- 3. `finite_field_traits<gf_ext_element<K, Irred>>` in lam-interop
  So NTT over extension fields works automatically when the degree is large 
  enough — relevant for Reed-Solomon and similar.

  The lambot-side dispatch bridge would look like `dispatch_gf_ext<N, F>` 
mirroring `dispatch_poly<N, F>` — parameterized on extension degree rather than 
polynomial degree
