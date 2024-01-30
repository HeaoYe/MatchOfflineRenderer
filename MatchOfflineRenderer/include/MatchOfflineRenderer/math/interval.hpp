#pragma once

#include <MatchOfflineRenderer/math/types.hpp>
#include <cmath>

namespace MatchOfflineRenderer::math {
    inline uint32_t float_to_bits(float rhs) {
        uint32_t dst;
        std::memcpy(&dst, &rhs, sizeof(uint32_t));
        return dst;
    }

    inline float bits_to_float(uint32_t rhs) {
        float dst;
        std::memcpy(&dst, &rhs, sizeof(float));
        return dst;
    }

    inline float next_float_up(float rhs) {
        if (std::isinf(rhs) && rhs > 0.0f) {
            return rhs;
        }
        if (rhs == -0.0f) {
            return 0.0f;
        }
        uint32_t ui = float_to_bits(rhs);
        if (rhs >= 0.0f) {
            ui ++;
        } else {
            ui --;
        }
        return bits_to_float(ui);
    }

    inline float next_float_down(float rhs) {
        if (std::isinf(rhs) && rhs < 0.0f)
            return rhs;
        if (rhs == 0.0f)
            rhs = -0.0f;
        uint32_t ui = float_to_bits(rhs);
        if (rhs > 0) {
            ui --;
        } else {
            ui ++;
        }
        return bits_to_float(ui);
    }

    inline Real add_round_up(Real lhs, Real rhs) {
        return next_float_up(lhs + rhs);
    }

    inline Real add_round_down(Real lhs, Real rhs) {
        return next_float_down(lhs + rhs);
    }

    inline Real sub_round_up(Real lhs, Real rhs) {
        return add_round_up(lhs, -rhs);
    }

    inline Real sub_round_down(Real lhs, Real rhs) {
        return add_round_down(lhs, -rhs);
    }

    inline Real mul_round_up(Real lhs, Real rhs) {
        return next_float_up(lhs * rhs);
    }

    inline Real mul_round_down(Real lhs, Real rhs) {
        return next_float_down(lhs * rhs);
    }

    inline Real div_round_up(Real lhs, Real rhs) {
        return next_float_up(lhs / rhs);
    }

    inline Real div_round_down(Real lhs, Real rhs) {
        return next_float_down(lhs / rhs);
    }

    inline Real sqrt_round_up(Real rhs) {
        return next_float_up(std::sqrt(rhs));
    }

    inline Real sqrt_round_down(Real rhs) {
        return next_float_down(std::sqrt(rhs));
    }

    inline Real fma_round_up(Real a, Real b, Real c) {
        return next_float_up(std::fma(a, b, c));
    }

    inline Real fma_round_down(Real a, Real b, Real c) {
        return next_float_down(std::fma(a, b, c));
    }

    struct Interval {
        Real low {}, high {};

        Interval() noexcept = default;

        Interval(Real low, Real high) noexcept : low(std::min(low, high)), high(std::max(low, high)) {}

        explicit Interval(Real rhs) noexcept : low(rhs), high(rhs) {}
        
        Real operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 2))
            if (i == 0) return low;
            return high;
        }

        Real upper_bound() const noexcept { return high; }

        Real lower_bound() const noexcept { return low; }

        Real mid_point() const noexcept { return (low + high) / 2; }

        Real width() const noexcept { return high - low; }

        
        explicit operator Real() const noexcept { return mid_point(); }

        bool exactly(Real rhs) const noexcept { return low == rhs && high == rhs; } 

        bool operator==(Real rhs) const noexcept { return exactly(rhs); }

        bool operator!=(Real rhs) const noexcept { return rhs < low || rhs > high; }

        bool operator==(const Interval &rhs) const noexcept {
            return low == rhs.low && high == rhs.high;
        }

        bool is_inside(Real rhs) const noexcept {
            return rhs >= lower_bound() && rhs <= upper_bound();
        }

        bool is_overlaps(const Interval &rhs) const noexcept {
            return lower_bound() <= rhs.upper_bound() &&
                upper_bound() >= rhs.lower_bound();
        }

        Interval operator-() const noexcept { return { -high, -low }; }
 
        Interval operator+(const Interval &rhs) const {
            return { add_round_down(low, rhs.low), add_round_up(high, rhs.high) };
        }
 
        Interval operator-(const Interval &rhs) const {
            return { sub_round_down(low, rhs.low), sub_round_up(high, rhs.high) };
        }

        Interval &operator+=(const Interval &rhs) {
            low = add_round_down(low, rhs.low);
            high = add_round_up(high, rhs.high);
            return *this;
        }

        Interval &operator-=(const Interval &rhs) {
            low = sub_round_down(low, rhs.low);
            high = sub_round_up(high, rhs.high);
            return *this;
        }

        Interval operator*(const Interval &rhs) const {
            Real lp[4] = { mul_round_down(low, rhs.low),  mul_round_down(high, rhs.low), mul_round_down(low, rhs.high), mul_round_down(high, rhs.high) };
            Real hp[4] = { mul_round_up(low, rhs.low),  mul_round_up(high, rhs.low), mul_round_up(low, rhs.high), mul_round_up(high, rhs.high) };
            return {
                std::min({ lp[0], lp[1], lp[2], lp[3] }),
                std::max({ hp[0], hp[1], hp[2], hp[3] })
            };
        }

        Interval operator/(const Interval &rhs) const {
            if (rhs.is_inside(0)) {
                return { -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity() };
            }
            Real lq[4] = { div_round_down(low, rhs.low), div_round_down(high, rhs.low), div_round_down(low, rhs.high), div_round_down(high, rhs.high) };
            Real hq[4] = { div_round_up(low, rhs.low), div_round_up(high, rhs.low), div_round_up(low, rhs.high), div_round_up(high, rhs.high) };
            return {
                std::min({ lq[0], lq[1], lq[2], lq[3] }),
                std::max({ hq[0], hq[1], hq[2], hq[3] })
            };
        }

        Interval &operator*=(const Interval &rhs) {
            *this = (*this) * rhs;
            return *this;
        }

        Interval &operator/=(const Interval &rhs) {
            *this = (*this) / rhs;
            return *this;
        }

        Interval &operator+=(Real rhs) { return *this += Interval(rhs); }

        Interval &operator-=(Real rhs) { return *this -= Interval(rhs); }

        Interval &operator*=(Real rhs) {
            if (rhs > 0) {
                *this = Interval(mul_round_down(rhs, low), mul_round_up(rhs, high));
            } else {
                *this = Interval(mul_round_down(rhs, high), mul_round_up(rhs, low));
            }
            return *this;
        }

        Interval &operator/=(Real rhs) {
            if (rhs > 0) {
                *this = Interval(div_round_down(low, rhs), div_round_up(high, rhs));
            } else {
                *this = Interval(div_round_down(high, rhs), div_round_up(low, rhs));
            }
            return *this;
        }

        Interval &square() {
            Real alow = std::abs(lower_bound()), ahigh = std::abs(upper_bound());
            if (alow > ahigh) {
                std::swap(alow, ahigh);
            }
            if (is_inside(0)) {
                low = 0;
                high = mul_round_up(ahigh, ahigh);
                return *this;
            }
            low = mul_round_down(alow, alow);
            high = mul_round_up(ahigh, ahigh);
            return *this;
        }

        inline static Interval from_value_and_error(Real value, Real error) {
            Interval result;
            if (error == 0) {
                result.low = value;
                result.high = value;
            } else {
                result.low = sub_round_down(value, error);
                result.high = add_round_up(value, error);
            }
            return result;
        }
    };
    
    inline Interval operator+(Real lhs, const Interval &rhs) {
        return Interval(lhs) + rhs;
    }
    
    inline Interval operator-(Real lhs, const Interval &rhs) {
        return Interval(lhs) - rhs;
    }

    inline Interval operator*(Real lhs, const Interval &rhs) {
        if (lhs > 0) {
            return { mul_round_down(lhs, rhs.lower_bound()), mul_round_up(lhs, rhs.upper_bound()) };
        }
        return { mul_round_down(lhs, rhs.upper_bound()), mul_round_up(lhs, rhs.lower_bound()) };
    }

    inline Interval operator/(Real lhs, const Interval &rhs) {
        if (rhs.is_inside(0)) {
            return { -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity() };
        }
        if (lhs > 0) {
            return { div_round_down(lhs, rhs.upper_bound()), div_round_up(lhs, rhs.lower_bound()) };
        }
        return { div_round_down(lhs, rhs.lower_bound()), div_round_up(lhs, rhs.upper_bound()) };
    }
    
    inline Interval operator+(const Interval &rhs, Real lhs) {
        return Interval(lhs) + rhs;
    }
    
    inline Interval operator-(const Interval &rhs, Real lhs) {
        return Interval(lhs) - rhs;
    }

    inline Interval operator*(const Interval &rhs, Real lhs) {
        if (lhs > 0) {
            return { mul_round_down(lhs, rhs.lower_bound()), mul_round_up(lhs, rhs.upper_bound()) };
        }
        return { mul_round_down(lhs, rhs.upper_bound()), mul_round_up(lhs, rhs.lower_bound()) };
    }

    inline Interval operator/(const Interval &rhs, Real lhs) {
        if (rhs.is_inside(0)) {
            return { -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity() };
        }
        if (lhs > 0) {
            return { div_round_down(lhs, rhs.upper_bound()), div_round_up(lhs, rhs.lower_bound()) };
        }
        return { div_round_down(lhs, rhs.lower_bound()), div_round_up(lhs, rhs.upper_bound()) };
    }

    inline Interval square(const Interval &rhs) {
        Real alow = std::abs(rhs.lower_bound()), ahigh = std::abs(rhs.upper_bound());
        if (alow > ahigh)
            std::swap(alow, ahigh);
        if (rhs.is_inside(0)) {
            return { 0, mul_round_up(ahigh, ahigh) };
        }
        return { mul_round_down(alow, alow), mul_round_up(ahigh, ahigh) };
    }

    template <>
    inline bool isnan<Interval>(Interval rhs) noexcept {
        return std::isnan(rhs.low) || std::isnan(rhs.high);
    }
}
