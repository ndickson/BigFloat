# BigFloat
An easy-to-use, C++17 class for high-precision floating-point arithmetic.

Here's an example computing sqrt(tau) (i.e. sqrt(2pi)) to 225 bits of precision, to ensure that it's very likely rounded to the closest possible double:

```
using namespace big_float;
BigFloat<7> f; // 7 means 7*32 + 1 = 225 mantissa bits
f.tau();       // Could alternatively copy variable constants::tau<7>
f.sqrt();
double result = (double)f;
```

Most of the functions can run in constexpr code, so can be evaluated at compile time, to avoid slow runtime performance.  The code isn't very optimized for runtime performance, but shouldn't be extremely slow.  The algorithms are intended for moderately high precsion, maybe up to 1024 bits, not millions of bits, so other libraries will be better in situations where very high precision is needed.

*NOTE:* It's not extensively tested, so please don't rely on the results for anything very important unless you verify results carefully.
