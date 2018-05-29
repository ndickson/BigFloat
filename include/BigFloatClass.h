#pragma once

// This file is for if the BigFloat class definition is needed,
// but the function implementations are not needed, e.g. if
// it's needed for a member variable declaration.

#include "BigFloatFwd.h" // This is just to ensure BigFloatFwd is consistent.
#include <stdint.h>
#include <stddef.h>

namespace big_float {

using int8 = int8_t;
using int16 = int16_t;
using int32 = int32_t;
using int64 = int64_t;
using uint8 = uint8_t;
using uint16 = uint16_t;
using uint32 = uint32_t;
using uint64 = uint64_t;

template<size_t THE_N>
struct BigFloat {
	static constexpr size_t N = THE_N;

	uint32 mantissa[N];
	int16 exponent;
	bool negative;

	static constexpr int16 EXP_INF_OR_NAN = 0x7FFF;
	static constexpr int16 EXP_ZERO = int16(0x8000);
	static constexpr uint32 MANTISSA_INF = 0;
	static constexpr uint32 MANTISSA_NAN = 1;

	BigFloat() noexcept = default;
	constexpr BigFloat(BigFloat&&) noexcept = default;
	constexpr BigFloat(const BigFloat&) noexcept = default;
	constexpr BigFloat& operator=(BigFloat&&) noexcept = default;
	constexpr BigFloat& operator=(const BigFloat&) noexcept = default;

	struct zero_init {
		constexpr zero_init() noexcept {};
		constexpr zero_init(const zero_init&) noexcept = default;
	};
	constexpr BigFloat(zero_init) : mantissa{0}, exponent{EXP_ZERO}, negative{false} {}
	struct pow2_init {
		constexpr pow2_init() noexcept {};
		constexpr pow2_init(const pow2_init&) noexcept = default;
	};
	constexpr BigFloat(int16 exponent_,pow2_init) noexcept : mantissa{0}, exponent{exponent_}, negative{false} {}
	struct tau_init {
		constexpr tau_init() noexcept {};
		constexpr tau_init(const tau_init&) noexcept = default;
	};
	constexpr BigFloat(int16 exponent_,tau_init) noexcept : mantissa{0}, exponent{0}, negative{false} {
		tau();
		*this <<= exponent_;
	}
	struct e_init {
		constexpr e_init() noexcept {};
		constexpr e_init(const e_init&) noexcept = default;
	};
	constexpr BigFloat(int16 exponent_,e_init) noexcept : mantissa{0}, exponent{0}, negative{false} {
		e();
		*this <<= exponent_;
	}

	constexpr explicit BigFloat(int64 v) noexcept;
	constexpr explicit BigFloat(uint64 v) noexcept;
	constexpr explicit BigFloat(int8 v) noexcept : BigFloat(int64(v)) {}
	constexpr explicit BigFloat(uint8 v) noexcept : BigFloat(uint64(v)) {}
	constexpr explicit BigFloat(int16 v) noexcept : BigFloat(int64(v)) {}
	constexpr explicit BigFloat(uint16 v) noexcept : BigFloat(uint64(v)) {}
	constexpr explicit BigFloat(int32 v) noexcept : BigFloat(int64(v)) {}
	constexpr explicit BigFloat(uint32 v) noexcept : BigFloat(uint64(v)) {}
	constexpr BigFloat& operator=(int64 v) noexcept;
	constexpr BigFloat& operator=(uint64 v) noexcept;
	constexpr BigFloat& operator=(int8 v) noexcept { *this = int64(v); return *this; }
	constexpr BigFloat& operator=(int16 v) noexcept { *this = int64(v); return *this; }
	constexpr BigFloat& operator=(int32 v) noexcept { *this = int64(v); return *this; }
	constexpr BigFloat& operator=(uint8 v) noexcept { *this = uint64(v); return *this; }
	constexpr BigFloat& operator=(uint16 v) noexcept { *this = uint64(v); return *this; }
	constexpr BigFloat& operator=(uint32 v) noexcept { *this = uint64(v); return *this; }
	/// NOTE: This can't be constexpr, because unfortunately,
	///       the C++ standard doesn't allow accessing the bits of a double.
	explicit BigFloat(double v) noexcept;
	/// Conversion from BigFloat with a different N.
	template<size_t OTHERN>
	constexpr explicit BigFloat(const BigFloat<OTHERN>& that) noexcept;
	/// Conversion from BigFloat with a different N.
	template<size_t OTHERN>
	constexpr BigFloat& operator=(const BigFloat<OTHERN>& that) noexcept;
	constexpr bool operator==(const BigFloat& other) const noexcept;
	constexpr bool operator!=(const BigFloat& other) const noexcept;
	constexpr bool operator<(const BigFloat& other) const noexcept;
	constexpr bool operator<=(const BigFloat& other) const noexcept;
	constexpr bool operator>(const BigFloat& other) const noexcept;
	constexpr bool operator>=(const BigFloat& other) const noexcept;
	constexpr bool isNaN() const noexcept;
	constexpr bool isInfinite() const noexcept;
	constexpr bool isFinite() const noexcept;
	constexpr bool isZero() const noexcept;
	/// NOTE: This can't be constexpr, because unfortunately,
	///       the C++ standard doesn't allow accessing the bits of a double.
	explicit operator double() const noexcept;
	/// NOTE: This can't be constexpr, because unfortunately,
	///       the C++ standard doesn't allow accessing the bits of a float.
	explicit operator float() const noexcept;
	constexpr void operator+=(const BigFloat& other) noexcept;
	constexpr BigFloat operator+(const BigFloat& other) const noexcept;
	/// Unary negation operator
	constexpr BigFloat operator-() const noexcept;
	/// Unary plus operator
	constexpr BigFloat operator+() const noexcept;
	constexpr void operator-=(const BigFloat& other) noexcept;
	constexpr BigFloat operator-(const BigFloat& other) const noexcept;

	/// Helper function used by operator+= and operator-=
	/// Computes sign(this)*(|this|+|other|),
	/// which is this+other if they have the same sign.
	constexpr void addSameSign(const BigFloat& other) noexcept;

	/// Helper function used by operator+= and operator-=
	/// Computes this-other, assuming this and other have the same sign as this.
	/// NOTE: This flips the sign of this if |this| < |other|,
	///       since that's what the subtraction would do.
	constexpr void subSameSign(const BigFloat& other) noexcept;

	/// Helper function used by subSameSign
	/// Computes this-other, assuming this and other have the same sign as this,
	/// and this is strictly larger than other, and neither is zero, infinite, or NaN.
	constexpr void subSameSignStage2(const BigFloat& other) noexcept;

	constexpr void operator*=(const BigFloat& other) noexcept;
	constexpr void operator*=(uint32 other) noexcept;
	constexpr void operator*=(int32 other) noexcept;
	constexpr BigFloat operator*(const BigFloat& other) const noexcept;
	constexpr BigFloat operator*(uint32 other) const noexcept;
	constexpr BigFloat operator*(int32 other) const noexcept;
	constexpr void operator/=(const BigFloat& other) noexcept;
	constexpr void operator/=(uint32 other) noexcept;
	constexpr void operator/=(int32 other) noexcept;
	constexpr BigFloat operator/(const BigFloat& other) const noexcept;
	constexpr BigFloat operator/(uint32 other) const noexcept;
	constexpr BigFloat operator/(int32 other) const noexcept;

	/// Divides by 2^bits
	constexpr void operator>>=(int16 bits) noexcept;
	constexpr BigFloat operator>>(int16 bits) const noexcept;
	/// Multiplies by 2^bits
	constexpr void operator<<=(int16 bits) noexcept;
	constexpr BigFloat operator<<(int16 bits) const noexcept;

	/// Assigns sqrt(other) to this
	constexpr void sqrt(const BigFloat& other) noexcept;

	/// Assigns tau (a.k.a. 2pi) to this
	constexpr void tau() noexcept;

	/// Assigns e (a.k.a. exp(1)) to this
	constexpr void e() noexcept;

	/// Assigns sin(x) to this
	/// WARNING: THIS IS NOT READY FOR USE, EXCEPT FOR SMALL x!!!
	constexpr void sin(const BigFloat& x) noexcept;
};

namespace constants {
/// Compile-time constant zero value, in case that's needed
template<size_t N>
static constexpr BigFloat<N> zero = BigFloat<N>(typename BigFloat<N>::zero_init());
/// Compile-time constant one value, in case that's needed
template<size_t N>
static constexpr BigFloat<N> one = BigFloat<N>(0,typename BigFloat<N>::pow2_init());
/// Compile-time constant negative one value, in case that's needed
template<size_t N>
static constexpr BigFloat<N> negative_one = -BigFloat<N>::one;
/// Compile-time constant tau, a.k.a. 2pi.
template<size_t N>
static constexpr BigFloat<N> tau = BigFloat<N>(0,typename BigFloat<N>::tau_init());
/// Compile-time constant pi, a.k.a. tau/2.
template<size_t N>
static constexpr BigFloat<N> pi = BigFloat<N>(-1,typename BigFloat<N>::tau_init());
/// Compile-time constant e, a.k.a. exp(1).
template<size_t N>
static constexpr BigFloat<N> e = BigFloat<N>(0,typename BigFloat<N>::e_init());
}

}
