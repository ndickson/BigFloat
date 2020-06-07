#pragma once

// This file is for if implementations of functions in the BigFloat class
// are needed.

#include "BigFloatClass.h"
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace big_float {

#if 0
struct DoubleInt {
	union {
		double d;
		int64 i;
	};
	DoubleInt() noexcept = default;
	constexpr DoubleInt(const DoubleInt&) noexcept = default;
	constexpr DoubleInt(DoubleInt&&) noexcept = default;
	constexpr DoubleInt& operator=(const DoubleInt&) noexcept = default;
	constexpr DoubleInt& operator=(DoubleInt&&) noexcept = default;
	constexpr DoubleInt(double d_) noexcept : d(d_) {}
	constexpr DoubleInt(int64 i_) noexcept : i(i_) {}
};
struct FloatInt {
	union {
		float f;
		int32 i;
	};
	FloatInt() noexcept = default;
	constexpr FloatInt(const FloatInt&) noexcept = default;
	constexpr FloatInt(FloatInt&&) noexcept = default;
	constexpr FloatInt& operator=(const FloatInt&) noexcept = default;
	constexpr FloatInt& operator=(FloatInt&&) noexcept = default;
	constexpr FloatInt(float f_) noexcept : f(f_) {}
	constexpr FloatInt(int32 i_) noexcept : i(i_) {}
};
#endif

constexpr static uint32 BitScanF64(uint64 v) {
	// NOTE: This doesn't use _BitScanForward64 (Visual Studio) or __builtin_ctzll (GCC/Clang)
	//       because _BitScanForward64 can't be called from constexpr code.
	if (v == 0) {
		return uint32(64);
	}
	bool has_32 = (uint32(v) != 0);
	uint32 count = uint32(!has_32)<<5;
	v = has_32 ? v : (v>>32);
	bool has_16 = ((v&0xFFFF) != 0);
	count += uint32(!has_16)<<4;
	v = has_16 ? v : (v>>16);
	bool has_8 = ((v&0xFF) != 0);
	count += uint32(!has_8)<<3;
	v = has_8 ? v : (v>>8);
	bool has_4 = ((v&0xF) != 0);
	count += uint32(!has_4)<<2;
	v = has_4 ? v : (v>>4);
	bool has_2 = ((v&0x3) != 0);
	count += uint32(!has_2)<<1;
	v = has_2 ? v : (v>>2);
	count += !(v&1);
	return count;
}
constexpr static uint32 BitScanR64(uint64 v) {
	// NOTE: This doesn't use _BitScanReverse64 (Visual Studio) or __builtin_clzll (GCC/Clang)
	//       because _BitScanReverse64 can't be called from constexpr code.
	if (v == 0) {
		return uint32(-1);
	}
	bool has_32 = ((v>>32) != 0);
	uint32 count = uint32(has_32)<<5;
	v = has_32 ? (v>>32) : (v);
	bool has_16 = ((v>>16) != 0);
	count += uint32(has_16)<<4;
	v = has_16 ? (v>>16) : (v);
	bool has_8 = ((v>>8) != 0);
	count += uint32(has_8)<<3;
	v = has_8 ? (v>>8) : v;
	bool has_4 = ((v>>4) != 0);
	count += uint32(has_4)<<2;
	v = has_4 ? (v>>4) : v;
	bool has_2 = ((v>>2) != 0);
	count += uint32(has_2)<<1;
	v = has_2 ? (v>>2) : v;
	count += uint32(v>>1);
	return count;
}
constexpr static uint32 BitScanR32(uint32 v) {
	// NOTE: This doesn't use _BitScanReverse (Visual Studio) or __builtin_clz (GCC/Clang)
	//       because _BitScanReverse can't be called from constexpr code.
	if (v == 0) {
		return uint32(-1);
	}
	bool has_16 = ((v>>16) != 0);
	uint32 count = uint32(has_16)<<4;
	v = has_16 ? (v>>16) : (v);
	bool has_8 = ((v>>8) != 0);
	count += uint32(has_8)<<3;
	v = has_8 ? (v>>8) : v;
	bool has_4 = ((v>>4) != 0);
	count += uint32(has_4)<<2;
	v = has_4 ? (v>>4) : v;
	bool has_2 = ((v>>2) != 0);
	count += uint32(has_2)<<1;
	v = has_2 ? (v>>2) : v;
	count += uint32(v>>1);
	return count;
}
#if 0
#ifdef __GNUC__
constexpr uint32 BitScanF64(uint64 v)
{
	if (v == 0) {
		return uint32(64);
	}
	return uint32(__builtin_ctzll((unsigned long long)v));
}
constexpr uint32 BitScanR32(uint32 v)
{
	if (v == 0) {
		return uint32(-1);
	}
	return uint32(31-__builtin_clz(v));
}
#elif defined(_MSC_VER)
constexpr uint32 BitScanF64(uint64 v)
{
	unsigned long top_bit = 64;
	if (!_BitScanForward64(&top_bit, v)) {
		return uint32(64);
	}
	return uint32(top_bit);
}
constexpr uint32 BitScanR32(uint32 v)
{
	unsigned long top_bit = 32;
	if (!_BitScanReverse(&top_bit, v)) {
		return uint32(-1);
	}
	return uint32(top_bit);
}
#else
#error "Unsupported compiler!  Sorry.  Please implement an intrinsic for BSF or BSR instructions."
#endif
#endif

template<size_t THE_N>
constexpr BigFloat<THE_N>::BigFloat(int64 v) noexcept : negative(false), exponent(0), mantissa{0} {
	*this = v;
}
template<size_t THE_N>
constexpr BigFloat<THE_N>::BigFloat(uint64 v) noexcept : negative(false), exponent(0), mantissa{0} {
	*this = v;
}

template<size_t THE_N>
constexpr BigFloat<THE_N>& BigFloat<THE_N>::operator=(int64 v) noexcept {
	negative = false;
	for (size_t i = 0; i < N; ++i) {
		mantissa[i] = 0;
	}
	if (v == 0) {
		exponent = EXP_ZERO;
		return *this;
	}
	if (v < 0) {
		negative = true;
		v = -v;
	}
	uint32 top_bit = BitScanR64(uint64(v));
	exponent = int16(top_bit);
	if (top_bit == 0) {
		return *this;
	}
	static_assert(N > 0, "N must be at least 1");
	if (top_bit <= 32) {
		if (N > 1) {
			mantissa[N-2] = 0;
		}
		mantissa[N-1] = uint32(v << (32-top_bit));
	}
	else {
		if (N > 1) {
			mantissa[N-2] = uint32(v << (64-top_bit));
		}
		mantissa[N-1] = uint32(uint64(v) >> (top_bit-32));
		if (N == 1 && ((uint64(v) >> (top_bit-33))&1)) {
			if ((++mantissa[N-1]) == 0) {
				++exponent;
			}
		}
	}
	return *this;
}

template<size_t THE_N>
constexpr BigFloat<THE_N>& BigFloat<THE_N>::operator=(uint64 v) noexcept {
	negative = false;
	for (size_t i = 0; i < N; ++i) {
		mantissa[i] = 0;
	}
	if (v == 0) {
		exponent = EXP_ZERO;
		return *this;
	}
	uint32 top_bit = BitScanR64(uint64(v));
	exponent = int16(top_bit);
	if (top_bit == 0) {
		return *this;
	}
	static_assert(N > 0, "N must be at least 1");
	if (top_bit <= 32) {
		if (N > 1) {
			mantissa[N-2] = 0;
		}
		mantissa[N-1] = uint32(v << (32-top_bit));
	}
	else {
		if (N > 1) {
			mantissa[N-2] = uint32(v << (64-top_bit));
		}
		mantissa[N-1] = uint32(v >> (top_bit-32));
		if (N == 1 && ((v >> (top_bit-33))&1)) {
			if ((++mantissa[N-1]) == 0) {
				++exponent;
			}
		}
	}
	return *this;
}

template<size_t THE_N>
BigFloat<THE_N>::BigFloat(double v) noexcept : negative(false), exponent(0), mantissa{0} {
	static_assert(sizeof(int64) == sizeof(double), "This code requires that int64 is the same size as double.");
	int64& i = *reinterpret_cast<int64*>(&v);
	negative = (i < 0);
	exponent = int16(((i >> 52) & 0x7FF) - 0x3FF);
	uint64 v_mantissa = uint64(i&((int64(1)<<52)-1));
	if (exponent == 0x400) {
		// Infinity or NaN
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = v_mantissa ? MANTISSA_NAN : MANTISSA_INF;
		for (size_t index = 1; index < N; ++index) {
			mantissa[index] = 0;
		}
		return;
	}
	else if (exponent == -0x3FF) {
		if (v_mantissa == 0) {
			// Zero
			exponent = EXP_ZERO;
			for (size_t index = 0; index < N; ++index) {
				mantissa[index] = 0;
			}
			return;
		}

		// Denormal: 0.mantissa x 2^(-0x3FE)
		uint32 top_bit = BitScanR64(v_mantissa);
		for (size_t index = 0; index < N-2; ++index) {
			mantissa[index] = 0;
		}
		if (top_bit <= 32) {
			if (N > 1) {
				mantissa[N-2] = 0;
			}
			mantissa[N-1] = uint32(v_mantissa << (32-top_bit));
		}
		else {
			if (N > 1) {
				mantissa[N-2] = uint32(v_mantissa << (64-top_bit));
			}
			mantissa[N-1] = uint32(v_mantissa >> (top_bit-32));
			if (N == 1 && ((v_mantissa >> (top_bit-33))&1)) {
				if ((++mantissa[N-1]) == 0) {
					// This increases the exponent by 1.
					// It won't overflow to infinity.
					++top_bit;
				}
			}
		}
		exponent = int16(-0x3FE - 52) + int16(top_bit);
	}
	if (N > 1) {
		for (size_t index = 0; index < N-2; ++index) {
			mantissa[index] = 0;
		}
		mantissa[N-2] = uint32(i << 12);
		mantissa[N-1] = uint32(i >> 20);
	}
	else {
		static_assert(N > 0, "N cannot be zero!");
		mantissa[0] = uint32(i >> 20);
		if (i&0x80000) {
			++mantissa[0];
			if (mantissa[0] == 0) {
				// NOTE: This won't overfloat to infinity, because we have 16 bits of exponent.
				++exponent;
			}
		}
	}
}
/// Conversion from BigFloat with a different N.
template<size_t THE_N>
template<size_t OTHERN>
constexpr BigFloat<THE_N>::BigFloat(const BigFloat<OTHERN>& that) noexcept : negative(that.negative), exponent(that.exponent), mantissa{0} {
	*this = that;
}
/// Conversion from BigFloat with a different N.
template<size_t THE_N>
template<size_t OTHERN>
constexpr BigFloat<THE_N>& BigFloat<THE_N>::operator=(const BigFloat<OTHERN>& that) noexcept {
	negative = that.negative;
	exponent = that.exponent;
	if (that.exponent == EXP_ZERO) {
		return *this;
	}
	if (that.exponent == EXP_INF_OR_NAN) {
		mantissa[0] = that.mantissa[0];
		return *this;
	}
	if (N >= OTHERN) {
		// Zero-initialize any low bits.
		for (size_t i = 0; i < N-OTHERN; ++i) {
			mantissa[i] = 0;
		}
		// More space (or equal), so just copy.
		for (size_t i = 0; i < OTHERN; ++i) {
			mantissa[i+(N-OTHERN)] = that.mantissa[i];
		}
		return *this;
	}
	// Less space, so copy and round.
	uint32 carry = (that.mantissa[(OTHERN-N)-1]>>31);
	for (size_t i = 0; i < N; ++i) {
		uint32 result = that.mantissa[i+(OTHERN-N)] + carry;
		mantissa[i] = result;
		carry &= (result == 0);
	}
	if (carry) {
		// If we carried all the way up, mantissa is zero,
		// and we need to go to the next power of two.
		// NOTE: This will automatically go to infinity if exponent was
		//       one less than EXP_INF_OR_NAN, because MANTISSA_INF is zero.
		++exponent;
	}
	return *this;
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator==(const BigFloat& other) const noexcept {
	if (exponent != other.exponent) {
		return false;
	}
	// Exponents are now equal

	if (exponent == EXP_ZERO) {
		// Both zero, even if different signs
		return true;
	}
	if (negative != other.negative) {
		return false;
	}
	// Now non-zero with same signs

	if (exponent == EXP_INF_OR_NAN) {
		// NaNs are never equal, but infinities can be equal
		return (mantissa[0] == MANTISSA_INF);
	}

	// Finite with same sign and exponent: compare mantissa
	for (size_t i = 0; i < N; ++i) {
		if (mantissa[i] != other.mantissa[i]) {
			return false;
		}
	}
	return true;
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator!=(const BigFloat& other) const noexcept {
	return !(*this == other);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator<(const BigFloat& other) const noexcept {
	if (isNaN() || other.isNaN()) {
		return false;
	}
	if (negative != other.negative) {
		return negative && (exponent != EXP_ZERO || other.exponent != EXP_ZERO);
	}
	// Same sign now

	if (exponent != other.exponent) {
		return negative != (exponent < other.exponent);
	}
	// Exponents are now equal

	if (exponent == EXP_INF_OR_NAN) {
		// Both infinity, and infinity is not less than infinity (though it is considered equal)
		return false;
	}
	if (exponent == EXP_ZERO) {
		return false;
	}
	// Now finite, non-zero, same sign, same exponent

	// Finite with same sign and exponent: compare mantissa
	for (size_t i = N-1; i < N; --i) {
		if (mantissa[i] != other.mantissa[i]) {
			return (mantissa[i] < other.mantissa[i]);
		}
	}
	// Equal values, not less
	return false;
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator<=(const BigFloat& other) const noexcept {
	if (isNaN() || other.isNaN()) {
		return false;
	}
	if (negative != other.negative) {
		// -zero is <= +zero, because it's considered equal
		// +zero is <= -zero, because it's considered equal
		return negative || (exponent == EXP_ZERO && other.exponent == EXP_ZERO);
	}
	// Same sign now

	if (exponent != other.exponent) {
		return negative != (exponent < other.exponent);
	}
	// Exponents are now equal

	if (exponent == EXP_INF_OR_NAN) {
		// Both infinity
		return true;
	}
	if (exponent == EXP_ZERO) {
		return true;
	}
	// Now finite, non-zero, same sign, same exponent

	// Finite with same sign and exponent: compare mantissa
	for (size_t i = N-1; i < N; --i) {
		if (mantissa[i] != other.mantissa[i]) {
			return (mantissa[i] < other.mantissa[i]);
		}
	}
	// Equal values
	return true;
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator>(const BigFloat& other) const noexcept {
	return !isNaN() && !other.isNaN() && !(*this <= other);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::operator>=(const BigFloat& other) const noexcept {
	return !isNaN() && !other.isNaN() && !(*this < other);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::isNaN() const noexcept {
	return (exponent == EXP_INF_OR_NAN) && (mantissa[0] != MANTISSA_INF);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::isInfinite() const noexcept {
	return (exponent == EXP_INF_OR_NAN) && (mantissa[0] == MANTISSA_INF);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::isFinite() const noexcept {
	return (exponent != EXP_INF_OR_NAN);
}
template<size_t THE_N>
constexpr bool BigFloat<THE_N>::isZero() const noexcept {
	return (exponent == EXP_ZERO);
}
template<size_t THE_N>
BigFloat<THE_N>::operator double() const noexcept {
	double output_double;
	static_assert(sizeof(int64) == sizeof(double), "This code requires that int64 is the same size as double.");
	int64& output_int = *reinterpret_cast<int64*>(&output_double);
	output_int = int64(negative)<<63;
	if (isNaN()) {
		// Quiet NaN
		output_int |= 0x7FF8000000000000ULL;
		return output_double;
	}
	int64 mantissa_bits = (int64(mantissa[N-1])<<20) | (mantissa[N-2]>>12);
	// Round up if bit after is 1
	// FIXME: Round half to even, (although it's unlikely to occur with full-precision random source data if N is >= 4).
	mantissa_bits += ((mantissa[N-2]>>11)&1);
	if (exponent >= 0x400 || (exponent == 0x3FF && (mantissa_bits==0x0010000000000000ULL))) {
		// Infinity
		output_int |= 0x7FF0000000000000ULL;
		return output_double;
	}
	if (exponent < -0x3FF || (exponent == -0x3FF && (mantissa_bits!=0x0010000000000000ULL))) {
		// Denormal; we currently round to zero
		// FIXME: Handle denormal doubles correctly.
		return output_double;
	}
	output_int |= int64(exponent+0x3FF)<<52;
	// NOTE: Adding mantissa_bits may add 1 to the exponent if the carry bit
	//       bubbled all the way up.
	output_int += mantissa_bits;
	return output_double;
}
template<size_t THE_N>
BigFloat<THE_N>::operator float() const noexcept {
	float output_float;
	static_assert(sizeof(int32) == sizeof(float), "This code requires that int32 is the same size as float.");
	int32& output_int = *reinterpret_cast<int32*>(&output_float);
	output_int = int32(negative)<<31;
	if (isNaN()) {
		// Quiet NaN
		output_int |= 0x7FC00000;
		return output_float;
	}
	int32 mantissa_bits = int32(mantissa[N-1]>>11);
	// Round up if bit after is 1
	// FIXME: Round half to even, (although it's unlikely to occur with full-precision random source data if N is >= 3).
	mantissa_bits += ((mantissa[N-1]>>10)&1);
	if (exponent >= 0x80 || (exponent == 0x7F && (mantissa_bits==0x00800000))) {
		// Infinity
		output_int |= 0x7F800000;
		return output_float;
	}
	if (exponent < -0x7F || (exponent == -0x7F && (mantissa_bits!=0x00800000))) {
		// Denormal; we currently round to zero
		// FIXME: Handle denormal floats correctly.
		return output_float;
	}
	output_int |= int32(exponent+0x7F)<<23;
	// NOTE: Adding mantissa_bits may add 1 to the exponent if the carry bit
	//       bubbled all the way up.
	output_int += mantissa_bits;
	return output_float;
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator+=(const BigFloat& other) noexcept {
	if (this == &other) {
		if (exponent != EXP_INF_OR_NAN && exponent != EXP_ZERO) {
			// finite_self+finite_self = 2*finite_self
			++exponent;
			if (exponent == EXP_INF_OR_NAN) {
				// Overflow to infinity, not NaN.
				mantissa[0] = MANTISSA_INF;
			}
		}
		return;
	}
	if (negative == other.negative) {
		addSameSign(other);
	}
	else {
		// +a + -b --> +a - +b
		// -a + +b --> -a - -b
		subSameSign(other);
	}
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator+(const BigFloat& other) const noexcept {
	BigFloat<THE_N> v(*this);
	v += other;
	return v;
}
/// Unary negation operator
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator-() const noexcept {
	BigFloat<THE_N> f(*this);
	f.negative = !f.negative;
	return f;
}
/// Unary plus operator
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator+() const noexcept {
	return *this;
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator-=(const BigFloat& other) noexcept {
	if (negative == other.negative) {
		if (this == &other) {
			if (exponent == EXP_INF_OR_NAN) {
				// NaN - NaN == NaN
				// infinity - infinity == NaN
				mantissa[0] = MANTISSA_NAN;
				return;
			}
			// Otherwise zero with same sign.
			// This should be equivalent to what subSameSign
			// would do with two identical BigFloats.
			exponent = EXP_ZERO;
			return;
		}
		subSameSign(other);
	}
	else {
		// +a - -b --> +a + +b
		// -a - +b --> -a + -b
		addSameSign(other);
	}
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator-(const BigFloat& other) const noexcept {
	BigFloat<THE_N> v(*this);
	v -= other;
	return v;
}
/// Helper function used by operator+= and operator-=
/// Computes sign(this)*(|this|+|other|),
/// which is this+other if they have the same sign.
template<size_t THE_N>
constexpr void BigFloat<THE_N>::addSameSign(const BigFloat& other) noexcept {
	if (other.exponent == EXP_ZERO) {
		// Adding zero does nothing.
		return;
	}
	if (exponent == EXP_ZERO) {
		// NOTE: Use this sign, since we're pretending that other has the same sign as this.
		bool old_negative = negative;
		*this = other;
		negative = old_negative;
		return;
	}
	if (exponent == EXP_INF_OR_NAN || other.exponent == EXP_INF_OR_NAN) {
		if (isNaN()) {
			// NaN + anything == NaN
			return;
		}
		if (other.isNaN()) {
			// anything + NaN == NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// infinity + anything_positive == infinity
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}

	uint32 new_mantissa[N+1] = {0}; // NOTE: Value is not used; here for constexpr only!
	new_mantissa[N] = 1;
	uint32 exponent_diff = 0; // NOTE: Value is not used; here for constexpr only!
	const uint32* source = nullptr; // NOTE: Value is not used; here for constexpr only!
	int16 max_exponent = 0; // NOTE: Value is not used; here for constexpr only!
	if (exponent < other.exponent) {
		if (int32(other.exponent) >= int32(exponent) + int32(32*N+2)) {
			// this is so small that it can't contribute to other.
			// NOTE: Use this sign, since we're pretending that other has the same sign as this.
			bool old_negative = negative;
			*this = other;
			negative = old_negative;
			return;
		}
		max_exponent = other.exponent;
		exponent_diff = uint32(other.exponent-exponent);
		source = mantissa;
		for (size_t i = 0; i < N; ++i) {
			new_mantissa[i] = other.mantissa[i];
		}
	}
	else {
		if (int32(exponent) >= int32(other.exponent) + int32(32*N+2)) {
			// other is so small that it can't contribute.
			return;
		}
		max_exponent = exponent;
		exponent_diff = uint32(exponent-other.exponent);
		source = other.mantissa;
		for (size_t i = 0; i < N; ++i) {
			new_mantissa[i] = mantissa[i];
		}
	}

	uint32 carry = 0; // NOTE: Value must be initialized for constexpr; used by one branch below.
	uint32 nonzero_below_carry = 0;
	if (exponent_diff == 0) {
		//carry = 0;
	}
	else if (exponent_diff == 32*N+1) {
		carry = 1;

		// Whether whole expclit part of mantissa is zero
		// might affect whether the code below needs to round to even.
		for (size_t i = 0; i < N; ++i) {
			if (source[i] != 0) {
				nonzero_below_carry = 1;
				break;
			}
		}
	}
	else {
		uint32 exponent_diff_m1 = (exponent_diff-1);
		size_t index(exponent_diff_m1>>5);
		size_t bit(exponent_diff_m1 & 0x1F);
		carry = (source[index]>>bit)&1;
		// Whether all lower bits of mantissa are zero
		// might affect whether the code below needs to round to even.
		nonzero_below_carry = ((source[index] & ~(uint32(-int32(1))<<bit)) != 0);
		if (!nonzero_below_carry) {
			for (size_t i = 0; i < index; ++i) {
				if (source[i] != 0) {
					nonzero_below_carry = 1;
					break;
				}
			}
		}
	}

	const bool round_half_to_even = carry && !nonzero_below_carry;
	const uint32 orig_carry = carry;

	size_t source_index0 = (exponent_diff>>5);
	size_t source_bit_shift = (exponent_diff & 0x1F);
	if (source_bit_shift != 0) {
		size_t i = 0;
		for (; source_index0 < N-1; ++i, ++source_index0) {
			uint64 value = uint64(new_mantissa[i])
				+ uint64((source[source_index0]>>source_bit_shift) | (source[source_index0+1]<<(32-source_bit_shift)))
				+ carry;
			carry = uint32(value>>32);
			new_mantissa[i] = uint32(value);
		}
		// Add in the implicit 1 from source.
		uint64 value = 0; // NOTE: Value is not used; here for constexpr only!
		if (source_index0 < N) {
			value = uint64(new_mantissa[i])
				+ uint64((source[source_index0]>>source_bit_shift) | (uint32(1)<<(32-source_bit_shift)))
				+ carry;
		}
		else {
			// NOTE: This case should only happen when exponent_diff == 32*N+1.
			//       We need this codepath just to avoid the out-of-bounds array access above.
			value = uint64(new_mantissa[i]) + carry;
		}
		carry = uint32(value>>32);
		new_mantissa[i] = uint32(value);

		// Propagate the carry up if needed
		for (++i; carry && i < N+1; ++i) {
			uint64 value = uint64(new_mantissa[i]) + carry;
			carry = uint32(value>>32);
			new_mantissa[i] = uint32(value);
		}
	}
	else {
		// No bit shift
		size_t i = 0;
		for (; source_index0 < N; ++i, ++source_index0) {
			uint64 value = uint64(new_mantissa[i])
				+ uint64(source[source_index0])
				+ carry;
			carry = uint32(value>>32);
			new_mantissa[i] = uint32(value);
		}
		// Add the implicit 1 from source into carry.
		++carry;

		for (; carry && i < N+1; ++i) {
			uint64 value = uint64(new_mantissa[i]) + carry;
			carry = uint32(value>>32);
			new_mantissa[i] = uint32(value);
		}
	}

	if (new_mantissa[N] == 1) {
		exponent = max_exponent;

		if (round_half_to_even) {
			// Just clear bit 0 of new_mantissa if rounding half to even,
			// since we previously just rounded half up.
			new_mantissa[0] &= ~uint32(1);
		}

		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = new_mantissa[i];
		}
	}
	else if (max_exponent+1 == EXP_INF_OR_NAN) {
		// Number became infinite
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
	}
	else {
		// new_mantissa[N] is 2 or 3.
		exponent = max_exponent+1;

		// Handle possibility of needing to round to even,
		// taking into account that if orig_carry is 1,
		// 1 was already added to new_mantissa for that carry.
		uint32 low_result_bits = (new_mantissa[0]&3);
		uint32 carry = 0;
		// If result_bits is 0, it would round to the current
		// value regardless of orig_carry or nonzero_below_carry.
		if (low_result_bits != 0) {
			// Remove the carry that was already added to new_mantissa.
			low_result_bits -= orig_carry;
			// Thinking of a 4-bit number consisting of 2 bits in low_result_bits,
			// 1 bit in orig_carry, and 1 bit in nonzero_below_carry,
			// 0000-0011 should round down to 0000 (new carry 0)
			// 0100 should round even down to 0000 (new carry 0)
			// 0101-0111 or should round up to 1000 (new carry 1)
			// 1000-1011 should round down to 1000 (new carry 0)
			// 1100 should round even up to 10000 (new carry 1)
			// 1101 should round up to 10000 (new carry 1)
			// 1110-1111 already rounded up, so aren't in this branch.
			carry = (low_result_bits == 1 && (orig_carry || nonzero_below_carry)) || (low_result_bits == 3);
		}
		for (size_t i = 0; i < N; ++i) {
			uint64 value = uint64((new_mantissa[i]>>1) | (new_mantissa[i+1]<<31)) + carry;
			carry = uint32(value>>32);
			mantissa[i] = uint32(value);
		}
		// Carry propagated all the way up, so 11.11111111... rounded up to 100.00000000...
		if (carry) {
			++exponent;
			if (exponent == EXP_INF_OR_NAN) {
				mantissa[0] = MANTISSA_INF;
			}
		}
	}
}
/// Helper function used by operator+= and operator-=
/// Computes this-other, assuming this and other have the same sign as this.
/// NOTE: This flips the sign of this if |this| < |other|,
///       since that's what the subtraction would do.
template<size_t THE_N>
constexpr void BigFloat<THE_N>::subSameSign(const BigFloat& other) noexcept {
	if (other.exponent == EXP_ZERO) {
		// Subtracting zero does nothing.
		return;
	}
	if (exponent == EXP_ZERO) {
		// Subtracting from zero negates the sign.
		// NOTE: Use this sign, since we're pretending that other has the same sign as this.
		bool new_negative = !negative;
		*this = other;
		negative = new_negative;
		return;
	}
	if (exponent == EXP_INF_OR_NAN || other.exponent == EXP_INF_OR_NAN) {
		if (isNaN()) {
			// NaN - anything == NaN
			return;
		}
		if (other.isNaN()) {
			// anything - NaN == NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		if (exponent == EXP_INF_OR_NAN) {
			if (other.exponent == EXP_INF_OR_NAN) {
				// infinity - infinity == NaN
				mantissa[0] = MANTISSA_NAN;
				return;
			}
			// infinity - anything_positive == infinity
			return;
		}
		// anything_positive - infinity == -infinity
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		negative = !negative;
		return;
	}

	bool is_this_smaller = (exponent < other.exponent);
	if (!is_this_smaller && exponent == other.exponent) {
		for (size_t i = N-1; true; --i) {
			if (mantissa[i] < other.mantissa[i]) {
				is_this_smaller = true;
				break;
			}
			else if (mantissa[i] > other.mantissa[i]) {
				break;
			}
			if (i == 0) {
				// Exactly equal, so subtraction is zero
				exponent = EXP_ZERO;
				return;
			}
		}
	}
	if (is_this_smaller) {
		BigFloat<THE_N> other2(other);
		other2.negative = !negative;
		other2.subSameSignStage2(*this);
		*this = other2;
	}
	else {
		subSameSignStage2(other);
	}
}
/// Helper function used by subSameSign
/// Computes this-other, assuming this and other have the same sign as this,
/// and this is strictly larger than other, and neither is zero, infinite, or NaN.
template<size_t THE_N>
constexpr void BigFloat<THE_N>::subSameSignStage2(const BigFloat& other) noexcept {
	if (int32(other.exponent) <= int32(exponent)-int32(32*N+3)) {
		// other is too small to affect this, so no change
		return;
	}
	uint32 exponent_diff = uint32(exponent - other.exponent);

	// We need up to 3 bits of other's mantissa below this mantissa to
	// round correctly, (i.e. round half to even),
	// where the lowest bit is 0 if everything below the others is zero,
	// and it is 1 if anything below the others is nonzero.
	// NOTE: Values must be initialized for constexpr; values used by some branches below.
	uint32 bitsBelow = 0;
	if (exponent_diff == 0) {
		//bitsBelow = 0;
	}
	else if (exponent_diff == 32*N+2) {
		uint32 combined = other.mantissa[0];
		for (size_t i = 1; i < N; ++i) {
			combined |= other.mantissa[i];
		}
		// 001c
		bitsBelow = uint32(2) + uint32(combined != 0);
	}
	else if (exponent_diff == 32*N+1) {
		constexpr uint32 noTopBitMask = uint32(0x7FFFFFFF);
		uint32 combined = other.mantissa[N-1] & noTopBitMask;
		for (size_t i = 0; i < N-1; ++i) {
			combined |= other.mantissa[i];
		}
		// 01xc
		bitsBelow = uint32(4) + uint32((other.mantissa[N-1] >> 31) << 1) + uint32(combined != 0);
	}
	else {
		// At least one bit of overlap between the mantissas,
		// and at least one bit of non-overlap.
		size_t index(exponent_diff>>5);
		size_t bit(exponent_diff & 0x1F);
		// Move one bit lower
		index -= (bit == 0);
		bit = (bit-1) & 0x1F;
		bitsBelow = ((other.mantissa[index]>>bit) & 1) << 2;
		if (exponent_diff > 1) {
			// Move one bit lower
			index -= (bit == 0);
			bit = (bit-1) & 0x1F;
			bitsBelow |= ((other.mantissa[index]>>bit) & 1) << 1;
			if (exponent_diff > 2) {
				// Check if any lower bits are nonzero
				index -= (bit == 0);
				bit = (bit-1) & 0x1F;
				uint32 topBitsMask = (uint32(2) << uint32(bit)) - 1;
				uint32 combined = (other.mantissa[index] & topBitsMask);
				for (size_t i = 0; i < index; ++i) {
					combined |= other.mantissa[i];
				}
				bitsBelow |= (combined != 0);
			}
		}
	}

	// As if bitsBelow (3 bits) had been subtracted from the corresponding zero in this mantissa.
	uint32 belowNewMantissa = uint32(-int32(bitsBelow)) & 7;
	// There's a carry (borrow) if anything nonzero was subtracted from zero.
	uint32 carry = uint32(bitsBelow != 0);

	uint32 new_mantissa[N+1] = {0}; // NOTE: Value is not used; here for constexpr only!
	for (size_t i = 0; i < N; ++i) {
		new_mantissa[i] = mantissa[i];
	}
	new_mantissa[N] = 1;

	// Do the main part of the subtraction.
	const uint32* source = other.mantissa;
	size_t source_index0 = (exponent_diff>>5);
	size_t source_bit_shift = (exponent_diff & 0x1F);
	if (source_bit_shift != 0) {
		size_t i = 0;
		for (; source_index0 < N-1; ++i, ++source_index0) {
			uint64 value = uint64(new_mantissa[i])
				- uint64((source[source_index0]>>source_bit_shift) | (source[source_index0+1]<<(32-source_bit_shift)))
				- carry;
			carry = uint32(value>>32)&1;
			new_mantissa[i] = uint32(value);
		}
		// Subtract in the implicit 1 from source.
		uint64 value = 0; // NOTE: Value is not used; here for constexpr only!
		if (source_index0 < N) {
			value = uint64(new_mantissa[i])
				- uint64((source[source_index0]>>source_bit_shift) | (uint32(1)<<(32-source_bit_shift)))
				- carry;
		}
		else {
			// NOTE: This case should only happen when exponent_diff == 32*N+1 or 32*N+2.
			//       We need this codepath just to avoid the out-of-bounds array access above.
			value = uint64(new_mantissa[i]) - carry;
		}
		carry = uint32(value>>32)&1;
		new_mantissa[i] = uint32(value);

		// Propagate the carry up if needed
		for (++i; carry && i < N+1; ++i) {
			uint64 value = uint64(new_mantissa[i]) - carry;
			carry = uint32(value>>32)&1;
			new_mantissa[i] = uint32(value);
		}
	}
	else {
		// No bit shift
		size_t i = 0;
		for (; source_index0 < N; ++i, ++source_index0) {
			uint64 value = uint64(new_mantissa[i])
				- uint64(source[source_index0])
				- carry;
			carry = uint32(value>>32)&1;
			new_mantissa[i] = uint32(value);
		}
		// Add the implicit 1 from source into carry.
		++carry;

		for (; carry && i < N+1; ++i) {
			uint64 value = uint64(new_mantissa[i]) - carry;
			carry = uint32(value>>32)&1;
			new_mantissa[i] = uint32(value);
		}
	}

	if (new_mantissa[N] == 1) {
		// Exponent stays the same.
		// Round up if applicable, (round half to even).
		uint32 odd = (new_mantissa[0] & 1);
		// Combine the lowest 2 bits, since only need to check
		// if nonzero below highest here.
		if ((belowNewMantissa & 4) != 0 && ((belowNewMantissa & 3) != 0 || odd)) {
			++new_mantissa[0];
			carry = (new_mantissa[0] == 0);
			// NOTE: This should never go out of bounds, or even increment the 1 at new_mantissa[N].
			for (size_t i = 1; carry; ++i) {
				++new_mantissa[i];
				carry = (new_mantissa[i] == 0);
			}
		}

		// Copy the mantissa.
		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = new_mantissa[i];
		}
	}
	else {
		// mantissa[N] should be 0, so exponent needs to decrease,
		// unless rounding undoes it.
		// NOTE: The only way for the exponent to be decreased by more than 1
		// is if exponent_diff is zero, in which case there are no bits below,
		// so we only need one extra bit, though it depends on the
		// 3 bitsBelowMantissa computed above.

		// First, round up, if applicable, (round half to even).
		bool odd = (belowNewMantissa & 4) != 0;
		uint32 bitBelow = 0;
		if ((belowNewMantissa & 2) != 0 && ((belowNewMantissa & 1) != 0 || odd)) {
			bitBelow = uint32(!odd);
			if (odd) {
				++new_mantissa[0];
				carry = (new_mantissa[0] == 0);
				// NOTE: This should never go out of bounds, but it may
				// increment the 0 at new_mantissa[N], making exactly a power of 2
				for (size_t i = 1; carry; ++i) {
					++new_mantissa[i];
					carry = (new_mantissa[i] == 0);
				}
				if (new_mantissa[N] == 1) {
					// Rounded back up to a power of 2 with the same exponent as before.
					for (size_t i = 0; i < N; ++i) {
						mantissa[i] = 0;
					}
					return;
				}
			}
		}

		// Find highest set bit.
		size_t i = N-1;
		for (; true; --i) {
			if (new_mantissa[i] != 0) {
				break;
			}
			if (i == 0) {
				// Subtracted to zero.  This can only happen if exponent_diff
				// is 0, so there are no bits below mantissa.
				exponent = EXP_ZERO;
				return;
			}
		}
		uint32 bit = BitScanR32(new_mantissa[i]);
		uint32 exponent_diff = uint32(32*(N-i)) - bit;
		int32 new_exponent = int32(exponent) - int32(exponent_diff);
		if (new_exponent <= int32(EXP_ZERO)) {
			// Underflow
			exponent = EXP_ZERO;
			return;
		}
		int16 old_exponent = exponent;
		exponent = new_exponent;
		size_t j = 0;
		if (bit != 0) {
			for (; j+1 <= i; ++j) {
				mantissa[N-1-j] = ((new_mantissa[i-j]<<(32-bit)) | (new_mantissa[i-j-1]>>bit));
			}
			// Add bitBelow bit in position corresponding with its original power,
			// since there's now room to represent it.
			mantissa[N-1-j] = ((new_mantissa[i-j]<<(32-bit)) | (bitBelow<<(32-bit-1)));
		}
		else {
			for (; j+1 <= i; ++j) {
				mantissa[N-1-j] = (new_mantissa[i-j-1]);
			}
			// Add bitBelow bit in position corresponding with its original power,
			// since there's now room to represent it.
			mantissa[N-1-j] = (bitBelow<<31);
		}
		for (++j; j < N; ++j) {
			mantissa[N-1-j] = 0;
		}
	}
}

template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator*=(const BigFloat& other) noexcept {
	// NOTE: Be careful not to write to any members of this until we're done reading
	//       any members of other, just in case they're the same object.

	negative ^= other.negative;

	// Handle existing infinities, NaNs, and zeros
	if (exponent == EXP_INF_OR_NAN) {
		// NaN times anything is NaN
		// Infinity times zero is NaN
		// Infinity times NaN is NaN
		if (mantissa[0] == MANTISSA_NAN || other.exponent == EXP_ZERO || (other.exponent == EXP_INF_OR_NAN && other.mantissa[0] == MANTISSA_NAN)) {
			// NaN
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// Infinity times something other than zero or NaN is infinity
		return;
	}
	else if (exponent == EXP_ZERO) {
		if (other.exponent == EXP_INF_OR_NAN) {
			// Zero times infinity or NaN is NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// Already zero
		return;
	}

	int32 new_exponent = exponent + other.exponent;

	constexpr size_t N2 = 2*N+1;
	uint32 new_mantissa[N2] = {0};
	// First, 1 times other.mantissa
	new_mantissa[N2-1] = 1;
	for (size_t i = 0; i < N; ++i) {
		new_mantissa[i+N] = other.mantissa[i];
	}
	//for (size_t i = 0; i < N; ++i) {
		//new_mantissa[i] = 0;
	//}
	for (size_t this_index = 0; this_index < N; ++this_index) {
		const uint64 this_value = uint64(mantissa[this_index]);
		uint32 carry = 0;
		for (size_t other_index = 0; other_index < N; ++other_index) {
			//  (2^32 - 1)^2 + (2^32 - 1) + (2^32 - 1)
			// = 2^64 - 2(2^32) + 1 + 2(2^32) - 2
			// = 2^64 - 1, so this is safe, albeit just barely
			uint64 result = this_value*other.mantissa[other_index] + carry + new_mantissa[this_index + other_index];
			new_mantissa[this_index + other_index] = uint32(result);
			carry = uint32(result>>32);
		}
		// Implicit 1 bit for other mantissa
		uint64 result = this_value*1 + carry + new_mantissa[this_index+N];
		new_mantissa[this_index+N] = uint32(result);
		carry = uint32(result>>32);
		// Make sure to propagate the carry up as far as necessary
		for (size_t other_index = N+1; carry; ++other_index) {
			uint64 result = carry + new_mantissa[this_index+other_index];
			new_mantissa[this_index+other_index] = uint32(result);
			carry = uint32(result>>32);
		}
	}
	if (new_mantissa[N2-1] > 1) {
		++new_exponent;
		uint32 carry = new_mantissa[N]&1;
		for (size_t i = 0; i < N; ++i) {
			uint64 result = uint64((new_mantissa[N+i]>>1) | (new_mantissa[N+i+1]<<31)) + carry;
			mantissa[i] = uint32(result);
			carry = uint32(result>>32);
		}
	}
	else {
		uint32 carry = (new_mantissa[N-1]>>31);
		for (size_t i = 0; i < N; ++i) {
			uint64 result = uint64(new_mantissa[N+i]) + carry;
			mantissa[i] = uint32(result);
			carry = uint32(result>>32);
		}
		if (carry) {
			// Carry bubbled all the way up.
			// Mantissa should be all zero right now,
			// though, so we don't need to shift it.
			++new_exponent;
		}
	}

	if (new_exponent >= EXP_INF_OR_NAN) {
		// Infinite
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = 0;
		return;
	}
	else if (new_exponent <= EXP_ZERO) {
		// Round to zero
		exponent = EXP_ZERO;
		return;
	}
	exponent = new_exponent;
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator*=(uint32 other) noexcept {
	if (other == 1) {
		// No change.
		return;
	}

	// Handle existing infinities, NaNs, and zeros
	if (exponent == EXP_INF_OR_NAN) {
		// NaN times anything is NaN
		// Infinity times zero is NaN
		if (mantissa[0] == MANTISSA_NAN || other == 0) {
			// NaN
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// Infinity times something greater than zero is infinity
		return;
	}
	else if (exponent == EXP_ZERO) {
		// Already zero
		return;
	}
	else if (other == 0) {
		exponent = EXP_ZERO;
		return;
	}

	// other is at least 2 at this point

	constexpr size_t N2 = N+2;
	uint32 new_mantissa[N2] = {0};
	for (size_t i = 0; i < N; ++i) {
		uint64 local_product = uint64(mantissa[i]) * other;
		uint32 local_product_low = uint32(local_product);
		uint32 local_product_high = uint32(local_product>>32);
		new_mantissa[i] += local_product_low;
		uint32 carry = (new_mantissa[i] < local_product_low);
		// NOTE: This addition should never produce a new carry, so we don't need uint64.
		new_mantissa[i+1] = local_product_high + carry;
	}
	// Handle implicit 1 bit at the top of the mantissa
	new_mantissa[N] += other;
	new_mantissa[N+1] = (new_mantissa[N] < other);

	if (new_mantissa[N+1]) {
		// Exactly 32 bits higher.
		int32 new_exponent = int32(exponent) + int32(32);
		if (new_exponent >= int32(EXP_INF_OR_NAN)) {
			// Overflow to infinity
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_INF;
			return;
		}
		uint32 carry = (new_mantissa[0]>>31);
		for (size_t i = 0; i < N; ++i) {
			uint64 result = uint64(new_mantissa[i+1]) + uint64(carry);
			mantissa[i] = uint32(result);
			carry = uint32(result>>32);
		}
		// NOTE: In this case of 32 bits higher, carry should never be 1 here.
#if 0
		if (carry) {
			// Rounded all the way up to new power of two.  mantissa is all zero.
			++new_exponent;
			if (new_exponent >= int32(EXP_INF_OR_NAN)) {
				// Overflow to infinity
				exponent = EXP_INF_OR_NAN;
				mantissa[0] = MANTISSA_INF;
				return;
			}
		}
#endif
		exponent = int16(new_exponent);
		return;
	}

	// NOTE: top_bit will be at least 1, because other is at least 2.
	uint32 top_bit = BitScanR32(new_mantissa[N]);
	int32 new_exponent = int32(exponent) + int32(top_bit);
	uint32 carry = (mantissa[0]>>(top_bit-1))&1;
	for (size_t i = 0; i < N; ++i) {
		uint64 result = uint64(uint32(new_mantissa[i]>>top_bit) | uint32(new_mantissa[i+1]<<(32-top_bit))) + uint64(carry);
		mantissa[i] = uint32(result);
		carry = uint32(result>>32);
	}
	if (carry) {
		// Rounded all the way up to new power of two.  mantissa is all zero.
		++new_exponent;
	}
	if (new_exponent >= int32(EXP_INF_OR_NAN)) {
		// Overflow to infinity
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}
	exponent = int16(new_exponent);
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator*=(int32 other) noexcept {
	if (other < 0) {
		negative = !negative;
		other = -other;
	}
	*this *= uint32(other);
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator*(const BigFloat& other) const noexcept {
	BigFloat<THE_N> v(*this);
	v *= other;
	return v;
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator*(uint32 other) const noexcept {
	BigFloat<THE_N> v(*this);
	v *= other;
	return v;
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator*(int32 other) const noexcept {
	BigFloat<THE_N> v(*this);
	v *= other;
	return v;
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator/=(const BigFloat& other) noexcept {
	negative ^= other.negative;

	if (exponent == EXP_ZERO) {
		if (other.exponent == EXP_ZERO || (other.exponent == EXP_INF_OR_NAN && other.mantissa[0] == MANTISSA_NAN)) {
			// 0/0 is NaN, 0/NaN is NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// 0/infinity is 0 and 0/finite is 0
		return;
	}
	if (exponent == EXP_INF_OR_NAN) {
		if (mantissa[0] == MANTISSA_NAN || other.exponent == EXP_INF_OR_NAN) {
			// NaN/anything is NaN, infinity/infinity is NaN, infinity/NaN is NaN
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// infinity/finite is infinity
		return;
	}
	if (other.exponent == EXP_ZERO) {
		// finite/0 is infinity
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}
	if (other.exponent == EXP_INF_OR_NAN) {
		if (other.mantissa[0] == MANTISSA_NAN) {
			// finite/NaN is NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// finite/infinity is 0
		exponent = EXP_ZERO;
		return;
	}
	if (this == &other) {
		// finite_self/finite_self == 1.0
		// NOTE: negative is already false from above.
		exponent = 0;
		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = 0;
		}
		return;
	}

	int32 new_exponent = int32(exponent) - int32(other.exponent);
	// new_exponent won't get any larger, so if it's EXP_ZERO or less, round to zero.
	if (new_exponent <= EXP_ZERO) {
		exponent = EXP_ZERO;
		return;
	}
	// new_exponent may get one smaller, so if it's EXP_INF_OR_NAN+1 or larger, round to infinity.
	// We have to check this later without the +1 if new_exponent is not decreased.
	if (new_exponent >= EXP_INF_OR_NAN+1) {
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}

	// Strictly greater than the top 33 bits of other's mantissa, (including the implicit 1),
	// but close enough that it eliminates as many bits as possible on each pass.
	const uint64 local_denominator = (0x100000000ULL | uint64(other.mantissa[N-1])) + 1;

	// This needs to be longer than N to improve the chances of having the correct rounding.
	// Having it be at least twice N is advisable, so that at least one full pass can be done to full accuracy
	constexpr uint32 M = 3*N+3;
	uint32 remaining[M] = {0};
	uint32 quotient[M] = {0};

	//for (size_t i = 0; i < M-N-1; ++i) {
		//remaining[i] = 0;
	//}
	for (size_t i = 0; i < N; ++i) {
		remaining[M-N-1+i] = mantissa[i];
	}
	remaining[M-1] = 1;

	size_t top_index = M-1;
	size_t top_bit = 0;

	// Repeat until remaining is small enough
	while (true) {
		// Take the top 64 bits of remaining.
		const uint64 local_numerator = (uint64(remaining[top_index])<<(63-top_bit)) | (uint64(remaining[top_index-1])<<(31-top_bit)) | ((top_bit!=31) ? (uint64(remaining[top_index-2])>>(top_bit+1)) : 0);
		// NOTE: The local_quotient should always fit in 32 bits, since denominator is at least 2^32+1, and numerator is strictly less than 2^64.
		//       Specifically, it's always in the range from 0x40000000 to 0xFFFFFFFF
		uint32 local_quotient = uint32(local_numerator/local_denominator);

		// Add the new contribution to the quotient
		uint32 quotient_high = (local_quotient>>(31-top_bit));
		uint32 carry = 0;
		if (top_bit != 31) {
			uint32 quotient_low = (local_quotient<<(top_bit+1));
			uint64 addition_low = uint64(quotient[top_index-1]) + quotient_low;
			quotient[top_index-1] = uint32(addition_low);
			carry = uint32(addition_low>>32);
		}
		uint64 addition = uint64(quotient[top_index]) + uint64(quotient_high) + carry;
		quotient[top_index] = uint32(addition);
		carry = uint32(addition>>32);
		for (size_t i = top_index+1; carry; ++i) {
			addition = uint64(quotient[i]) + carry;
			quotient[i] = uint32(addition);
			carry = uint32(addition>>32);
		}

		// Subtract local_quotient times the full denominator from remaining.
		for (size_t i = 0; i < N; ++i) {
			const uint64 local_product = uint64(local_quotient)*uint64(other.mantissa[i]);
			// For the implicit 1 that would be at i==N, local_product would be local_quotient,
			// in which case, bit 31 would correspond with top_index and top_bit.
			// Since i is less than N, here, bit 31 corresponds with top_index-N+i and top_bit.
			const uint32 local_product_low = (top_bit!=31) ? uint32(local_product<<(top_bit+1)) : uint32(0);
			const uint32 local_product_mid = uint32(local_product>>(31-top_bit));
			const uint32 local_product_high = uint32(local_product>>(63-top_bit));

			uint32 carry = 0;
			if (top_bit != 31) {
				uint64 sub_low = uint64(remaining[top_index-N+i-1]) - uint64(local_product_low);
				remaining[top_index-N+i-1] = uint32(sub_low);
				carry = uint32((sub_low>>32)&1);
			}
			uint64 sub_mid = uint64(remaining[top_index-N+i]) - uint64(local_product_mid) - carry;
			remaining[top_index-N+i] = uint32(sub_mid);
			carry = uint32((sub_mid>>32)&1);
			if (local_product_high || carry) {
				uint64 sub_high = uint64(remaining[top_index-N+i+1]) - uint64(local_product_high) - carry;
				remaining[top_index-N+i+1] = uint32(sub_high);
				carry = uint32((sub_high>>32)&1);
				for (size_t j = top_index-N+i+2; carry; ++j) {
					uint64 sub = uint64(remaining[j]) - carry;
					remaining[j] = uint32(sub);
					carry = uint32((sub>>32)&1);
				}
			}
		}
		{
			// Remember to subtract local_quotient times implicit 1 from remaining,
			// to finish the above subtraction
			const uint64 local_product = uint64(local_quotient)*uint64(1);
			const uint32 local_product_low = (top_bit!=31) ? uint32(local_product<<(top_bit+1)) : uint32(0);
			const uint32 local_product_mid = uint32(local_product>>(31-top_bit));
			// NOTE: local_product_high would be zero here.
			uint32 carry = 0;
			if (top_bit != 31) {
				uint64 sub_low = uint64(remaining[top_index-1]) - uint64(local_product_low);
				remaining[top_index-1] = uint32(sub_low);
				carry = uint32((sub_low>>32)&1);
			}
			uint64 sub_mid = uint64(remaining[top_index]) - uint64(local_product_mid) - carry;
			remaining[top_index] = uint32(sub_mid);
			// NOTE: top_index was the highest index with at least one bit set,
			//       so we don't need to go higher.
		}

		// Find the new top_index and top_bit
		while (top_index > N && remaining[top_index]==0) {
			--top_index;
		}
		if (top_index == N) {
			// Code above accesses remaining[top_index-N-1], so top_index must be at least N+1,
			// else we have to stop.
			break;
		}

		top_bit = BitScanR32(remaining[top_index]);
	}

	size_t quotient_carry_index = M-1-N-1;

	// The high index of quotient should either be 0 or 1.
	// If it's 0, we need to reduce new_exponent.
	if (quotient[M-1] == 0) {
		// Round quotient up if bit 30 of element M-1-N-1 is 1
		uint32 carry = uint32(quotient[quotient_carry_index]>>30)&1;

		// Bits almost line up, so copy to mantissa with extra bit included from below.
		for (size_t i = 0; i < N; ++i) {
			uint64 value = uint64((quotient[M-1-N+i]<<1) | (quotient[M-1-N+i-1]>>31)) + carry;
			mantissa[i] = uint32(value);
			carry = uint32(value>>32);
		}
		if (carry) {
			// Rounded up to power of two, so we no longer have to decrease new_exponent.
			exponent = int16(new_exponent);
			return;
		}

		--new_exponent;
		if (new_exponent <= EXP_ZERO) {
			// Round to zero
			exponent = EXP_ZERO;
			return;
		}
		exponent = int16(new_exponent);
	}
	else {
		// Round quotient up if bit 31 of element M-1-N-1 is 1
		uint32 carry = uint32(quotient[quotient_carry_index]>>31);
		while (carry) {
			++quotient_carry_index;
			uint32 v = quotient[quotient_carry_index] + 1;
			quotient[quotient_carry_index] = v;
			carry = (v==0);
		}

		exponent = int16(new_exponent);

		// Bits line up, so copy to mantissa
		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = quotient[M-1-N+i];
		}
	}
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator/=(uint32 denominator) noexcept {
	if (denominator == 1) {
		// finite/1 is unchanged, infinity/1 is infinity, NaN/1 is NaN
		return;
	}
	if (exponent == EXP_ZERO) {
		if (denominator == 0) {
			// 0/0 is NaN
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		// 0/finite is 0
		return;
	}
	if (exponent == EXP_INF_OR_NAN) {
		// NaN/anything is NaN
		// infinity/finite is infinity and infinity/0 is infinity
		return;
	}
	if (denominator == 0) {
		// finite/0 is infinity
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}

	// other is now at least 2, so exponent will strictly decrease.

	// Unlike in the case of dividing by BigFloat, we don't have to add 1 to denominator here,
	// since local_denominator is exactly the denominator.

	// This needs to be longer than N to improve the chances of having the correct rounding.
	// We don't need quite as much precision as the case of dividing by BigFloat.
	constexpr uint32 M = 2*N+3;
	uint32 remaining[M] = {0};

	//for (size_t i = 0; i < M-N-1; ++i) {
	//remaining[i] = 0;
	//}
	for (size_t i = 0; i < N; ++i) {
		remaining[M-N-1+i] = mantissa[i];
	}
	remaining[M-1] = 1;

	// Because denominator is at least 2, and remaining[M-1] is 1, so quotient can be one index shorter.
	// NOTE: Unlike in the case of dividing by BigFloat, quotient here is treated more like an integer
	//       in the same power-of-two scale as remaining.
	uint32 quotient[M-1] = {0};

	size_t top_index = M-1;
	uint32 top_bit = 0;

	// Repeat until remaining is small enough
	while (true) {
		// Take the top 64 bits of remaining.
		const uint64 local_numerator = (uint64(remaining[top_index])<<(63-top_bit)) | (uint64(remaining[top_index-1])<<(31-top_bit)) | ((top_bit!=31) ? (uint64(remaining[top_index-2])>>(top_bit+1)) : 0);

		// NOTE: Unlike in the case of dividing by BigFloat, this can be larger than 32 bits, since we didn't shift the denominator.
		const uint64 local_quotient = local_numerator/denominator;

		// Add the new contribution to the quotient, remembering to account for the shift we didn't do.
		const uint32 local_quotient_low = (top_bit!=31) ? uint32(local_quotient<<(top_bit+1)) : uint32(0);
		const uint32 local_quotient_mid = uint32(local_quotient>>(31-top_bit));
		const uint32 local_quotient_high = uint32(local_quotient>>(63-top_bit));
		uint32 carry = 0;
		if (local_quotient_low) {
			uint64 result = uint64(quotient[top_index-2]) + uint64(local_quotient_low);
			quotient[top_index-2] = uint32(result);
			carry = uint32(result>>32);
		}
		if (local_quotient_mid || carry) {
			uint64 result = uint64(quotient[top_index-1]) + uint64(local_quotient_mid) + uint64(carry);
			quotient[top_index-1] = uint32(result);
			carry = uint32(result>>32);
		}
		if (local_quotient_high || carry) {
			uint64 result = uint64(quotient[top_index]) + uint64(local_quotient_high) + uint64(carry);
			quotient[top_index] = uint32(result);
			carry = uint32(result>>32);
		}
		for (size_t i = top_index+1; carry; ++i) {
			uint64 result = uint64(quotient[i]) + uint64(carry);
			quotient[i] = uint32(result);
			carry = uint32(result>>32);
		}

		// Subtract local_quotient times denominator from remaining.
		const uint64 local_product = local_quotient*denominator;
		const uint32 local_product_low = (top_bit!=31) ? uint32(local_product<<(top_bit+1)) : uint32(0);
		const uint32 local_product_mid = uint32(local_product>>(31-top_bit));
		const uint32 local_product_high = uint32(local_product>>(63-top_bit));
		carry = 0;
		if (local_product_low) {
			uint64 result = uint64(remaining[top_index-2]) - uint64(local_product_low);
			remaining[top_index-2] = uint32(result);
			carry = uint32(result>>32)&1;
		}
		if (local_product_mid || carry) {
			uint64 result = uint64(remaining[top_index-1]) - uint64(local_product_mid) - uint64(carry);
			remaining[top_index-1] = uint32(result);
			carry = uint32(result>>32)&1;
		}
		if (local_product_high || carry) {
			uint64 result = uint64(remaining[top_index]) - uint64(local_product_high) - uint64(carry);
			remaining[top_index] = uint32(result);
			// NOTE: carry should always be zero here, since remaining should always
			//       decrease by at least several bits each time, so we don't need to go up higher.
			//carry = uint32(result>>32)&1;
		}

		// Find the new top_index and top_bit
		while (top_index > 2 && remaining[top_index]==0) {
			--top_index;
		}
		if (top_index == 2) {
			// Code above accesses remaining[top_index-N-1], so top_index must be at least N+1,
			// else we have to stop.
			break;
		}

		top_bit = BitScanR32(remaining[top_index]);
	}

	// 2^32 / (2^32 - 1) quotient is 1, and 2^32 / 2 quotient is 2^31,
	// so only quotient[M-2] matters for finding top bit.
	top_bit = BitScanR32(quotient[M-2]);

	// Always lose at least 1 bit, since denominator is at least 2.
	int32 new_exponent = int32(exponent) - int32(32-top_bit);

	// Copy mantissa, rounding up for any carry bit
	if (top_bit == 0) {
		uint32 carry = (quotient[M-3-N]>>31);
		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = quotient[M-2-N+i] + carry;
			carry &= (mantissa[i] == 0);
		}
		if (carry) {
			// mantissa is all zero; we rounded up to next power of two;
			// still can't be as high as original, though.
			++new_exponent;
		}
	}
	else {
		uint32 carry = (quotient[M-2-N]>>(top_bit-1))&1;
		for (size_t i = 0; i < N; ++i) {
			mantissa[i] = ((quotient[M-2-N+i]>>top_bit) | (quotient[M-2-N+i+1]<<(32-top_bit))) + carry;
			carry &= (mantissa[i] == 0);
		}
		if (carry) {
			// mantissa is all zero; we rounded up to next power of two;
			// still can't be as high as original, though.
			++new_exponent;
		}
	}

	if (new_exponent <= EXP_ZERO) {
		exponent = EXP_ZERO;
		return;
	}
	exponent = int16(new_exponent);
}
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator/=(int32 denominator) noexcept {
	if (denominator < 0) {
		negative = !negative;
		denominator = -denominator;
	}
	*this /= uint32(denominator);
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator/(const BigFloat& other) const noexcept {
	BigFloat<THE_N> v(*this);
	v /= other;
	return v;
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator/(uint32 other) const noexcept {
	BigFloat<THE_N> v(*this);
	v /= other;
	return v;
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator/(int32 other) const noexcept {
	BigFloat<THE_N> v(*this);
	v /= other;
	return v;
}

/// Divides by 2^bits
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator>>=(int16 bits) noexcept {
	if (exponent == EXP_ZERO || exponent == EXP_INF_OR_NAN) {
		// Zero times anything finite is zero.
		// Infinity times anything finite is infinity.
		// NaN time anything is NaN.
		return;
	}
	int32 new_exponent = int32(exponent)-int32(bits);
	if (new_exponent <= int32(EXP_ZERO)) {
		// Round to zero.
		exponent = EXP_ZERO;
		return;
	}
	if (new_exponent >= int32(EXP_INF_OR_NAN)) {
		// Round to infinity.
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}
	exponent = int16(new_exponent);
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator>>(int16 bits) const noexcept {
	BigFloat<THE_N> v(*this);
	v >>= bits;
	return v;
}
/// Multiplies by 2^bits
template<size_t THE_N>
constexpr void BigFloat<THE_N>::operator<<=(int16 bits) noexcept {
	if (exponent == EXP_ZERO || exponent == EXP_INF_OR_NAN) {
		// Zero times anything finite is zero.
		// Infinity times anything finite is infinity.
		// NaN time anything is NaN.
		return;
	}
	int32 new_exponent = int32(exponent)+int32(bits);
	if (new_exponent <= int32(EXP_ZERO)) {
		// Round to zero.
		exponent = EXP_ZERO;
		return;
	}
	if (new_exponent >= int32(EXP_INF_OR_NAN)) {
		// Round to infinity.
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}
	exponent = int16(new_exponent);
}
template<size_t THE_N>
constexpr BigFloat<THE_N> BigFloat<THE_N>::operator<<(int16 bits) const noexcept {
	BigFloat<THE_N> v(*this);
	v <<= bits;
	return v;
}

template<size_t THE_N>
constexpr void BigFloat<THE_N>::roundToNumBits(size_t num_mantissa_bits) noexcept {
	// If it's intrinsically rounded enough, there's nothing to do.
	constexpr size_t max_mantissa_bits = N*DATA_TYPE_BITS + 1;
	if (num_mantissa_bits >= max_mantissa_bits) {
		return;
	}

	// Rounding infinity, NaN, or zero, does nothing.
	if (exponent == EXP_INF_OR_NAN || exponent == EXP_ZERO) {
		return;
	}

	// It doesn't make sense to round to zero bits of mantissa, since there's
	// always the implicit 1.
	if (num_mantissa_bits <= 0) {
		num_mantissa_bits = 1;
	}

	size_t rounding_bit_index = (max_mantissa_bits-1 - num_mantissa_bits);
	const size_t rounding_index = rounding_bit_index >> LOG_DATA_TYPE_BITS;
	rounding_bit_index &= (DATA_TYPE_BITS-1);
	const size_t kept_index = rounding_index + (rounding_bit_index == DATA_TYPE_BITS-1);
	const size_t kept_bit_index = (rounding_bit_index + 1) & (DATA_TYPE_BITS-1);
	const DataType rounding_bit = (mantissa[rounding_index]>>rounding_bit_index) & 1;
	const DataType kept_bit = (kept_index == N) ? 1 : ((mantissa[kept_index]>>kept_bit_index) & 1);

	// Zero out low bits and check if any below the rounding bit were nonzero.
	bool nonzero_below_rounding = false;
	for (size_t index = 0; index < rounding_index; ++index) {
		nonzero_below_rounding |= (mantissa[index] != 0);
		mantissa[index] = 0;
	}
	const DataType zero_mask = (~DataType(0)) >> (DATA_TYPE_BITS-1-rounding_bit_index);
	const DataType below_rounding_mask = (zero_mask >> 1);
	nonzero_below_rounding |= ((mantissa[rounding_index] & below_rounding_mask) != 0);
	mantissa[rounding_index] &= ~zero_mask;

	// Thinking of a 3-bit number consisting of the lowest "kept" bit,
	// the "rounding bit", and 1 bit for nonzero_below_rounding,
	// 000-001 should round down to 000 (new carry 0)
	// 010 should round even down to 000 (new carry 0)
	// 011 or should round up to 100 (new carry 1)
	// 100-101 should round down to 100 (new carry 0)
	// 110 should round even up to 1000 (new carry 1)
	// 111 should round up to 1000 (new carry 1)
	DataType carry = rounding_bit & (kept_bit | DataType(nonzero_below_rounding));

	if (carry == 0) {
		// Carry is zero, so nothing to do.
		return;
	}

	if (kept_index == N) {
		++exponent;
		if (exponent == EXP_INF_OR_NAN) {
			mantissa[0] = MANTISSA_INF;
		}
		return;
	}

	// Propagate the carry up
	carry <<= kept_bit_index;
	size_t index = kept_index;
	do {
		DataType value = mantissa[index];
		value += carry;
		carry = (value == 0);
		mantissa[index] = value;
		++index;
		if (index == N && carry != 0) {
			++exponent;
			if (exponent == EXP_INF_OR_NAN) {
				mantissa[0] = MANTISSA_INF;
			}
			return;
		}
	} while (carry != 0);
}

/// Assigns sqrt(other) to this
template<size_t THE_N>
constexpr void BigFloat<THE_N>::sqrt(const BigFloat& other) noexcept {
	if (other.exponent == EXP_ZERO) {
		// sqrt(+0) = +0
		// sqrt(-0) = -0
		negative = other.negative;
		exponent = EXP_ZERO;
		return;
	}
	if (other.negative) {
		// sqrt(-finite) = +NaN
		negative = false;
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_NAN;
		return;
	}

	constexpr size_t BIGN = (3*N+1)/2 + 4;
	BigFloat<BIGN> y(other);
	BigFloat<BIGN> xi(y);

	// Scale x to be between 1 and 2, and y to be between 1 and 4,
	// to avoid issues with temporary underflow during solve.
	xi.exponent = 0;
	y.exponent = (other.exponent&1);

	// We use mantissa shifted back 1 with the extra bit being the low exponent bit
	// as an approximation, since that is the result for linear interpolation between powers of two.
	// NOTE: Don't use the sqrt function in the standard library for an initial guess,
	//       just in case there are differences between platforms.
	xi.mantissa[BIGN-N-1] = (xi.mantissa[BIGN-N]<<31);
	for (size_t i = 0; i < N-1; ++i) {
		xi.mantissa[BIGN-N+i] = (xi.mantissa[BIGN-N+i]>>1) | (xi.mantissa[BIGN-N+i+1]<<31);
	}
	xi.mantissa[BIGN-1] = (xi.mantissa[BIGN-1]>>1) | (uint32(other.exponent&1)<<31);

	// Iterate until the difference on some iteration is small.
	// The difference on one iteration is always larger than the remaining true difference.
	int16 diff_exponent = 0;
	for (size_t iteration = 0; iteration < 20 && diff_exponent > -32*int16(BIGN-2) && diff_exponent > EXP_ZERO; ++iteration) {

		// x ~~ xi + (y-yi)/(dy/dx)
		//    = xi + (y - xi*xi)/(2*xi)
		BigFloat<BIGN> yi(xi);
		yi *= xi;
		BigFloat<BIGN> diff(y);
		diff -= yi;
		BigFloat<BIGN> deriv(xi);
		deriv <<= 1;
		diff /= deriv;
		diff_exponent = diff.exponent;

		xi += diff;
	}

	// The exponent will be floor(exponent/2), since
	// if 2^(2k) <= x < 2^(2(k+1)), that means that
	// 2^k <= sqrt(x) < 2^(k+1),
	// except possibly taking into account rounding, since the
	// slope is lower for larger values.
	xi <<= (other.exponent>>1);

	// NOTE: sqrt brings values closer to 1, so we shouldn't have to worry about
	//       the final value being zero or infinite.
	*this = BigFloat<N>(xi);
}

/// Assigns tau (a.k.a. 2pi) to this
template<size_t THE_N>
constexpr void BigFloat<THE_N>::tau() noexcept {
	constexpr uint32 known_mantissa_bits[12] = {
		0x921FB544,
		0x42D18469,
		0x898CC517,
		0x01B839A2,
		0x52049C11,
		0x14CF98E8,
		0x04177D4C,
		0x76273644,
		0xA29410F3,
		0x1C6809BB,
		0xDF2A3367,
		0x9A748636
		// NOTE: High bit of next number would be zero.
	};
	if (N <= 12) {
		negative = false;
		exponent = 2;
		for (size_t i = 0; i < N; ++i) {
			mantissa[N-1-i] = known_mantissa_bits[i];
		}
		// NOTE: High bit of next number would be zero for case of N==12, so we can skip this.
		if (N < 12 && (known_mantissa_bits[N]&0x80000000)) {
			// NOTE: None of the values in known_mantissa_bits are 0xFFFFFFFF, so there's no carry.
			++mantissa[0];
		}
		return;
	}

	// For higher N, we need to compute it, using Machin's formula:
	// tau/8 = 4*atan(1/5) - atan(1/239)
	// atan(x) = x - (x^3)/3 + (x^5)/5 - (x^7)/7 + (x^9)/9 - ...

	// log_5(2^32) = 13.7816..., and the denominators are increasing too,
	// so 14*N terms should be plenty for the atan(1/5).
	constexpr size_t nterms_5 = 14*N + 4;

	// Compute at higher precision to reduce roundoff error,
	// particularly because we'll be summing in forward order.
	constexpr size_t BIGN = (3*N)/2+5;

	BigFloat<BIGN> power_5(0,typename BigFloat<BIGN>::pow2_init());
	power_5 /= uint32(5);
	BigFloat<BIGN> sum(power_5);
	//BigFloat<BIGN> recip_25(power_5);
	//recip_25 *= power_5;
	for (size_t termi = 3; termi < nterms_5; termi += 2) {
		//power_5 *= recip_25;
		power_5 /= uint32(25);
		BigFloat<BIGN> term(power_5);
		//term /= BigFloat<BIGN>(uint64(termi));
		term /= uint32(termi);
		term.negative = ((termi>>1)&1);
		sum += term;
	}
	// Multiply by 4
	sum.exponent += 2;

	// log_239(2^32) = 4.05..., and the denominators are increasing too,
	// so (9*N)/2 terms should be plenty for the atan(1/239).
	constexpr size_t nterms_239 = (9*N)/2 + 4;
	BigFloat<BIGN> power_239(0,typename BigFloat<BIGN>::pow2_init());
	power_239 /= uint32(239);
	sum -= power_239;
	//BigFloat<BIGN> recip_239_2(power_239);
	//recip_239_2 *= power_239;
	for (size_t termi = 3; termi < nterms_239; termi += 2) {
		//power_239 *= recip_239_2;
		power_239 /= uint32(239*239);
		BigFloat<BIGN> term(power_239);
		//term /= BigFloat<BIGN>(uint64(termi));
		term /= uint32(termi);
		term.negative = ((termi>>1)&1);
		sum -= term;
	}

	// Multiply by 8, since Machin's formula computes tau/8.
	sum.exponent += 3;

	// Reduce precision back down
	*this = sum;
}

/// Assigns e (a.k.a. exp(1)) to this
template<size_t THE_N>
constexpr void BigFloat<THE_N>::e() noexcept {
	constexpr uint32 known_mantissa_bits[12] = {
		0x5BF0A8B1,
		0x45769535,
		0x5FB8AC40,
		0x4E7A79E3,
		0xB1738B07,
		0x9C5A6D2B,
		0x53C26C82,
		0x28C867F7,
		0x99273B9C,
		0x49367DF2,
		0xFA5FC6C6,
		0xC618EBB1
		// NOTE: High bit of next number would be 1.
	};
	if (N <= 12) {
		negative = false;
		exponent = 1;
		for (size_t i = 0; i < N; ++i) {
			mantissa[N-1-i] = known_mantissa_bits[i];
		}
		// NOTE: High bit of next number would be one for case of N==12, so it has a carry.
		if (N == 12 || (known_mantissa_bits[N]&0x80000000)) {
			// NOTE: None of the values in known_mantissa_bits are 0xFFFFFFFF, so there's no carry.
			++mantissa[0];
		}
		return;
	}

	// Compute e using a Pade approximant of exp(1/(2^m)),
	// raised to the power of 2^m via repeated squaring.
	// The Pade approximant we're using is f8(x)/f8(-x), where f8(x)
	// = 518918400
	// + 259459200*x
	// +  60540480*x^2
	// +   8648640*x^3
	// +    831600*x^4
	// +     55440*x^5
	// +      2520*x^6
	// +        72*x^7
	// +           x^8
	// Since x will be a power of 1/2, and we can choose the power to avoid certain
	// carry cases, we can easily construct each in a single pass though mantissa.

	constexpr size_t BIGN = 2*N+4;

	// x will implicitly be 1/(2^m).
	// Every increment of m gives us about 17 more bits of accuracy,
	// and an m of 0 will give a bit more than 48 bits of accuracy.
	// The +16 is to round up, just in case.
	// Choosing 32*BIGN bits is overkill, since we don't have enough
	// precision to get that accuracy, but plays it safe.
	// NOTE: To be able to construct numerator and denominator exactly,
	//       8*m+28 must be <= 32*BIGN.  That's always satisfied,
	//       since that's equivalent to m <= (32*BIGN-28)/8, and with the
	//       choice below, m <= (32*BIGN-32)/16, but if it's
	//       ever increased to be more overkill, be careful.
	constexpr size_t m = (32*BIGN - 48 + 16)/17;

	// The exponent doesn't actually matter, since numerator and denominator
	// will have the same exponent, as long as m is at least 4; (there's
	// an overflow in the numerator if m is 3 or less).
	BigFloat<BIGN> numerator(0,BigFloat<BIGN>::pow2_init());
	BigFloat<BIGN> denominator(0,BigFloat<BIGN>::pow2_init());
	numerator.mantissa[BIGN-1]   = 0xEEE11000; // This is (518918400ULL<<4) - (1ULL<<32)
	denominator.mantissa[BIGN-1] = 0xEEE11000; //
	size_t full_bit = 32*BIGN - 28 - m;
	constexpr uint32 coeffs[8] = {
		259459200, 60540480, 8648640, 831600, 55440, 2520, 72, 1
	};
	for (size_t i = 0; i < 8; ++i) {
		const size_t index = full_bit/32;
		const size_t bit = full_bit - 32*index;
		const uint32 coeff = coeffs[i];
		const uint32 coeff_low = (coeff << bit);
		const uint32 coeff_high = bit ? (coeff >> (32-bit)) : 0;

		// Add to numerator
		{
			size_t local_index = index;
			const uint32 addition_low = numerator.mantissa[local_index] + coeff_low;
			numerator.mantissa[local_index] = addition_low;
			uint32 carry = (addition_low < coeff_low);
			uint32 coeff_add_high = coeff_high + carry;
			while (coeff_add_high) {
				++local_index;
				uint32 addition = numerator.mantissa[local_index] + coeff_add_high;
				numerator.mantissa[local_index] = addition;
				coeff_add_high = (addition < coeff_add_high);
			}
		}

		// Update denominator
		if (i&1) {
			// Even power (odd iteration), so addition
			size_t local_index = index;
			const uint32 addition_low = denominator.mantissa[local_index] + coeff_low;
			denominator.mantissa[local_index] = addition_low;
			uint32 carry = (addition_low < coeff_low);
			uint32 coeff_add_high = coeff_high + carry;
			while (coeff_add_high) {
				++local_index;
				uint32 addition = denominator.mantissa[local_index] + coeff_add_high;
				denominator.mantissa[local_index] = addition;
				coeff_add_high = (addition < coeff_add_high);
			}
		}
		else {
			// Odd power (even iteration), so subtraction
			size_t local_index = index;
			const uint32 subtraction_low = denominator.mantissa[local_index] - coeff_low;
			denominator.mantissa[local_index] = subtraction_low;
			uint32 carry = (coeff_low != 0) && (subtraction_low >= uint32(-int32(coeff_low)));
			uint32 coeff_sub_high = coeff_high + carry;
			while (coeff_sub_high) {
				++local_index;
				uint32 subtraction = denominator.mantissa[local_index] - coeff_sub_high;
				denominator.mantissa[local_index] = subtraction;
				// NOTE: Don't need to check if coeff_sub_high is non-zero, since we already checked.
				coeff_sub_high = (subtraction >= uint32(-int32(coeff_sub_high)));
			}
		}

		full_bit -= m;
	}

	numerator /= denominator;

	// Raise approximation of e^(1/(2^m)) to the power of 2^m, to approximate e^1
	for (size_t i = 0; i < m; ++i) {
		numerator *= numerator;
	}
	// Reduce the precision back down
	*this = numerator;
}

/// Assigns sin(x) to this
/// WARNING: THIS IS NOT READY FOR USE, EXCEPT FOR SMALL x!!!
template<size_t THE_N>
constexpr void BigFloat<THE_N>::sin(const BigFloat& x) noexcept {
	BigFloat<THE_N> terms[100];
	terms[0] = x;
	BigFloat<THE_N> x2 = x;
	x2 *= x;
	for (size_t i = 1; i < 100; ++i) {
		terms[i] = terms[i-1];
		terms[i] *= x2;
		terms[i] /= uint32((2*i)*(2*i+1));
		terms[i].negative = !terms[i].negative;
	}
	// Sum from smallest to largest, to reduce roundoff error.
	*this = terms[99];
	for (size_t i = 98; i < 100; --i) {
		*this += terms[i];
	}
}

/// Assigns cos(x) to this
/// WARNING: THIS IS NOT READY FOR USE, EXCEPT FOR SMALL x!!!
template<size_t THE_N>
constexpr void BigFloat<THE_N>::cos(const BigFloat& x) noexcept {
	BigFloat<THE_N> terms[100];
	terms[0] = 1;
	BigFloat<THE_N> x2 = x;
	x2 *= x;
	for (size_t i = 1; i < 100; ++i) {
		terms[i] = terms[i-1];
		terms[i] *= x2;
		terms[i] /= uint32((2*i-1)*(2*i));
		terms[i].negative = !terms[i].negative;
	}
	// Sum from smallest to largest, to reduce roundoff error.
	*this = terms[99];
	for (size_t i = 98; i < 100; --i) {
		*this += terms[i];
	}
}

/// Assigns exp(x) to this
template<size_t THE_N>
constexpr void BigFloat<THE_N>::exp(const BigFloat& x) noexcept {
	if (x.exponent == EXP_INF_OR_NAN) {
		if (x.mantissa[0] == MANTISSA_INF && x.negative) {
			// e^-infinity = +0
			negative = false;
			exponent = EXP_ZERO;
			return;
		}

		// If +infinity, output +infinity.  If NaN, output NaN.
		*this = x;
		return;
	}

	constexpr size_t BIGN = THE_N*4 + 1;
	BigFloat<BIGN> reduced(x);
	size_t numSquares = 0;
	constexpr ExponentType maxExponent = -ExponentType(DATA_TYPE_BITS*N)/4;
	// Reduce the exponent (i.e. divide by a power of two) so that
	// x is close enough to zero that the truncated Taylor series will be close enough.
	if (reduced.exponent > maxExponent) {
		numSquares = size_t(reduced.exponent - maxExponent);
		reduced.exponent = maxExponent;
	}

	// sum = 1 + x + (x^2)/2 + (x^3)/6 + (x^4)/24
	BigFloat<BIGN> sum(0, BigFloat<BIGN>::pow2_init());
	BigFloat<BIGN> current(reduced);
	sum += current;
	current *= reduced;
	current >>= 1;
	sum += current;
	current *= reduced;
	current /= 3;
	sum += current;
	current *= reduced;
	current /= 4;
	sum += current;

	// Repeatedly square to reintroduce the removed power of two.
	for (size_t i = 0; i < numSquares; ++i) {
		sum *= sum;
	}

	*this = BigFloat<THE_N>(sum);
}

/// Assigns ln(x) to this
template<size_t THE_N>
constexpr void BigFloat<THE_N>::ln(const BigFloat& x) noexcept {
	if (x.isZero()) {
		// ln(0) is negative infinity
		negative = true;
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_INF;
		return;
	}
	// NOTE: This is checked after checking for zero, since -0 and +0 are treated the same.
	if (x.negative) {
		// ln(x) is NaN for negative x
		negative = false;
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_NAN;
		return;
	}
	if (x.exponent == EXP_INF_OR_NAN) {
		// If infinity, output infinity.  If NaN, output NaN.
		*this = x;
		return;
	}

	constexpr size_t BIGN = THE_N*4 + 1;
	BigFloat<BIGN> reduced(x);
	
	// Repeatedly square root x until it is close enough to 1.
	size_t numSquareRoots = 0;
	while (reduced.exponent < -1 || reduced.exponent > 0) {
		reduced.sqrt(reduced);
		++numSquareRoots;
	}
	if (reduced.exponent == 0) {
		// TODO: This might not be a good heuristic for "close enough to 1".
		while (reduced.mantissa[BIGN-1] != 0) {
			reduced.sqrt(reduced);
			++numSquareRoots;
		}
	}
	else { // reduced.exponent == -1
		// TODO: This might not be a good heuristic for "close enough to 1".
		while (reduced.mantissa[BIGN-1] != ~DataType(0)) {
			reduced.sqrt(reduced);
			++numSquareRoots;
		}
	}

	// Compute a power series approximation of ln(x) for x close to 1.
	// It should converge faster than the Taylor series for ln(x) around 1,
	// and unlike the Taylor series, converges for all x > 0.
	// It's based on that ln(x) = ln(1 + (x-1)/(x+1)) - ln(1 - (x-1)/(x+1))
	// and expanding, combining, and simplifying the Taylor series for each of those logarithms.
	const BigFloat<BIGN> factor((reduced - constants::one<BIGN>)/(reduced + constants::one<BIGN>));
	const BigFloat<BIGN> factor2 = factor*factor;
	BigFloat<BIGN> sum(constants::one<BIGN>);
	BigFloat<BIGN> currentPower(factor2);
	// TODO: Is this enough terms in all cases, given the heuristic for "close enough" above?
	constexpr size_t numTerms = THE_N + 1;
	for (size_t i = 0; i < numTerms; ++i) {
		sum += currentPower / uint32(2*i + 1);
		currentPower *= currentPower;
	}
	sum += currentPower / uint32(2*numTerms + 1);

	sum *= factor;
	sum <<= 1;

	// Multiply by 2*numSquareRoots.
	sum *= uint32(2*numSquareRoots);

	*this = sum;
}

template<size_t THE_N>
constexpr void BigFloat<THE_N>::pow(const BigFloat& base, const BigFloat& expon) noexcept {
	if (base.isNaN() || expon.isNaN()) {
		negative = false;
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_NAN;
		return;
	}
	if (base.isZero()) {
		if (expon.isZero()) {
			// 0^0 is NaN
			negative = false;
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = MANTISSA_NAN;
			return;
		}
		if (expon.exponent == EXP_INF_OR_NAN && expon.negative) {
			// (+0)^-infinity = (+infinity)^+infinity = +infinity
			// (-0)^-infinity = (-infinity)^+infinity = NaN, since the sign is indeterminate
			negative = false;
			exponent = EXP_INF_OR_NAN;
			mantissa[0] = expon.negative ? MANTISSA_NAN : MANTISSA_INF;
			return;
		}

		// 0^exponent is 0
		negative = base.negative;
		exponent = EXP_ZERO;
		return;
	}
	if (expon.isZero()) {
		// base^0 is 1
		negative = false;
		exponent = 0;
		for (size_t i = 0; i < THE_N; ++i) {
			mantissa[i] = 0;
		}
		return;
	}

	if (base.negative) {
		// Only valid if expon is an integer.
		// FIXME: Implement this!!!

		// Otherwise NaN
		negative = false;
		exponent = EXP_INF_OR_NAN;
		mantissa[0] = MANTISSA_NAN;
		return;
	}

	// Simple case
	constexpr size_t BIGN = (THE_N*3)/2 + 1;
	BigFloat<BIGN> lnBase;
	lnBase.ln(BigFloat<BIGN>(base));

	BigFloat<BIGN> result;
	result.exp(lnBase*BigFloat<BIGN>(expon));

	*this = BigFloat<THE_N>(result);
}

}
