#pragma once

#include <MatchOfflineRenderer/commons.hpp>

#include <cstdint>

// 基础数学类型, 基础计算函数, 精确计算函数
namespace MatchOfflineRenderer::math {
    using Int = int32_t;
    // 实数
    // using Real = float;
    using Real = double;

    struct CompensatedReal {
        CompensatedReal(Real v, Real err = 0) : v(v), err(err) {}
        explicit operator float() const { return v + err; }
        explicit operator double() const { return double(v) + double(err); }
        Real v {}, err {};
    };

    static constexpr Real pi = Real { 3.14159265358979323846 };
    static constexpr Real inv_pi = Real { 0.31830988618379067154 };
    static constexpr Real inv2_pi = Real { 0.15915494309189533577 };
    static constexpr Real inv4_pi = Real { 0.07957747154594766788 };
    static constexpr Real pi_div2 = Real { 1.57079632679489661923 };
    static constexpr Real pi_div4 = Real { 0.78539816339744830961 };
    static constexpr Real sqrt2 = Real { 1.41421356237309504880 };

    inline Real radians(Real deg) noexcept { return (pi / 180) * deg; }
    
    inline Real degrees(Real rad) noexcept { return (180 / pi) * rad; }

    inline Real safe_asin(Real x) noexcept { return std::asin(std::clamp<Real>(x, -1, 1)); }
    
    inline Real safe_acos(Real x) noexcept { return std::acos(std::clamp<Real>(x, -1, 1)); }

    inline Real safe_sqrt(Real x) noexcept { return std::sqrt(std::max<Real>(x, 0)); }

    inline Real lerp(Real lhs, Real rhs, Real t) noexcept { return (1 - t) * lhs + t * rhs; }

    template <typename Ta, typename Tb, typename Tc, typename Td>
    inline auto difference_of_products(Ta a, Tb b, Tc c, Td d) noexcept {
        auto cd = c * d;
        auto difference_of_products_ = std::fma(a, b, -cd);
        auto error = std::fma(-c, d, cd);
        return difference_of_products_ + error;
    }

    inline CompensatedReal precise_product(Real lhs, Real rhs) noexcept {
        return { lhs * rhs, std::fma(lhs, rhs, lhs * rhs) };
    }

    inline CompensatedReal precise_sum(Real lhs, Real rhs) noexcept {
        Real sum = lhs + rhs;
        Real delta = sum - lhs;
        return { sum, (lhs - (sum - delta)) + (rhs - delta) };
    }

    namespace internel {
        template <typename T>
        inline CompensatedReal precise_inner_product(T lhs, T rhs) noexcept {
            return precise_product(lhs, rhs);
        }

        template <typename T, typename ...Ts>
        inline CompensatedReal precise_inner_product(T lhs, T rhs, Ts... args) noexcept {
            CompensatedReal lhs_rhs = precise_product(lhs, rhs);
            CompensatedReal args_p = precise_inner_product(args...);
            CompensatedReal sum = precise_sum(lhs_rhs.v, args_p.v);
            return { sum.v, lhs_rhs.err + (args_p.err + sum.err) };
        }
    }

    template <typename... Ts>
    inline std::enable_if_t<std::conjunction_v<std::is_arithmetic<Ts>...>, Real> precise_inner_product(Ts... args) noexcept {
        CompensatedReal ip = internel::precise_inner_product(args...);
        return Real { ip };
    }

    inline uint16_t real_encode_uint16(Real rhs) noexcept {
        return std::round(std::clamp<Real>((rhs + 1) / 2, 0, 1) * Real { 65535 });
    }

    template <typename T>
    bool isnan(T t) noexcept {
        return std::isnan(t);
    }

    template <typename T, typename C>
    constexpr T evaluate_polynomial(T t, C c) { return c; }

    template <typename T, typename C, typename... Args>
    constexpr T evaluate_polynomial(T t, C c, Args... args) {
        return std::fma(t, evaluate_polynomial(t, args...), c);
    }

    inline Real sigmoid(Real rhs) {
        if (std::isinf(rhs)) {
            return rhs > 0 ? 1 : 0;
        }
        return 0.5 + rhs / (2 * std::sqrt(1 + rhs * rhs));
    };
}
