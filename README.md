# lam
LA Math: A budding, modularized math library in c++23, and beyond.

In early construction phase. That being said, `lam` should build with the right combination of `cmake` commands, and `ninja`. If for some reason it does not build for you, please contact me somehow!

Currently `lam` consists of a few modules:
- [`lam.bitvector`](https://www.github.com/colinrford/bitvector)
- [`lam.concepts`](https://www.github.com/colinrford/concepts)
- [`lam.ctbignum`](https://www.github.com/colinrford/ctbignum)
- [`lam.lebesgue`](https://www.github.com/colinrford/lebesgue)
- [`lam.linearalgebra`](https://www.github.com/colinrford/linearalgebra)
- [`lam.polynomial_nttp`](https://www.github.com/colinrford/polynomial_nttp)
- [`lam.symbols`](https://www.github.com/colinrford/symbols)

At this time there is a small amount of interoperability between `lam` modules:
- `lam.bitvector` comes with a `std::ranges::sized_range` constructor and an `export_words` function, which are compatible with both `bigint` and `ZqElement` from `lam.ctbignum`
- root-finding interoperability between `lam.ctbignum` and `lam.polynomial_nttp`, try at your own risk! for now :)
- `lam.lebesgue` uses `lam.linearalgebra` and `lam.polynomial_nttp` internally

At this time `lam` does not have a universal license, so seek out the individual licenses in each of the subprojects.
