#pragma once

#include <MatchOfflineRenderer/math/tuple.hpp>
#include <MatchOfflineRenderer/math/interval.hpp>

#include <sstream>
#include <fmt/core.h>
#include <fmt/format.h>


namespace MatchOfflineRenderer {
    // 2D向量和3D向量
    
    template <typename T>
    struct Point2;

    template <typename T>
    struct Vector2 : public Tuple2<Vector2, T> {
        using Tuple2<Vector2, T>::x;
        using Tuple2<Vector2, T>::y;

        Vector2() noexcept = default;
        
        Vector2(T x, T y) noexcept : Tuple2<Vector2, T>(x, y) {}

        template <typename R>
        explicit Vector2(const Point2<R> &rhs) noexcept;

        template <typename R>
        explicit Vector2(const Vector2<R> &rhs) noexcept : Tuple2<Vector2, T>(T { rhs.x }, T { rhs.y }) {}

        T length_squa() const noexcept {
            return x * x + y * y;
        }

        tuple_length_t<T> length() const noexcept {
            return std::sqrt<tuple_length_t<T>>(x * x + y * y);
        }

        Vector2 &normalize() noexcept {
            auto inv_length = tuple_length_t<T> { 1 } / length();
            x *= inv_length;
            y *= inv_length;
            return *this;
        }

        T dot(const Vector2 &rhs) const noexcept {
            return x * rhs.x + y * rhs.y;
        }

        T abs_dot(const Vector2 &rhs) const noexcept {
            return std::abs(x * rhs.x + y * rhs.y);
        }

        Real angle(const Vector2 &rhs) const noexcept {
            if (dot(rhs) < 0) {
                return pi - 2 * safe_asin(length(*this + rhs) / 2);
            } else {
                return 2 * safe_asin(length(rhs - *this) / 2);
            }
        }

        Vector2 &gram_schmidt(const Vector2 &rhs) noexcept {
            *this -= rhs * dot(rhs);
            return *this;
        }
    };

    template <typename T>
    inline T length_squa(const Vector2<T> &rhs) {
        return rhs.x * rhs.x + rhs.y * rhs.y;
    }

    template <typename T>
    inline tuple_length_t<T> length(const Vector2<T> &rhs) {
        return std::sqrt<tuple_length_t<T>>(rhs.x * rhs.x + rhs.y * rhs.y);
    }

    template <typename T>
    inline Vector2<T> normalize(const Vector2<T> &rhs) {
        auto inv_length = tuple_length_t<T> { 1 } / rhs.length();
        return { rhs.x * inv_length, rhs.y * inv_length };
    }

    template <typename T>
    inline T dot(const Vector2<T> &lhs, const Vector2<T> &rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y;
    }

    template <typename T>
    inline T abs_dot(const Vector2<T> &lhs, const Vector2<T> &rhs) {
        return std::abs(lhs.x * rhs.x + lhs.y * rhs.y);
    }

    template <typename T>
    inline T angle(const Vector2<T> &lhs, const Vector2<T> &rhs) {
        if (dot(lhs, rhs) < 0) {
            return pi - 2 * safe_asin(length(lhs + rhs) / 2);
        } else {
            return 2 * safe_asin(length(rhs - lhs) / 2);
        }
    }

    template <typename T>
    inline Vector2<T> gram_schmidt(const Vector2<T> &lhs, const Vector2<T> &rhs) {
        return lhs - rhs * dot(lhs, rhs);
    }

    template <typename T>
    struct Point3;
    template <typename T>
    struct Normal3;

    template <typename T>
    struct Vector3 : public Tuple3<Vector3, T> {
        using Tuple3<Vector3, T>::x;
        using Tuple3<Vector3, T>::y;
        using Tuple3<Vector3, T>::z;

        Vector3() noexcept = default;

        Vector3(T x, T y, T z) noexcept : Tuple3<Vector3, T>(x, y, z) {}

        template <typename R>
        explicit Vector3(const Point3<R> &rhs) noexcept;

        template <typename R>
        explicit Vector3(const Normal3<R> &rhs) noexcept;

        template <typename R>
        explicit Vector3(const Vector3<R> &rhs) noexcept : Tuple3<Vector3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

        T length_squa() const noexcept {
            return x * x + y * y + z * z;
        }

        tuple_length_t<T> length() const noexcept {
            return std::sqrt(x * x + y * y + z * z);
        }

        Vector3 &normalize() noexcept {
            auto inv_length = tuple_length_t<T> { 1 } / length();
            x *= inv_length;
            y *= inv_length;
            z *= inv_length;
            return *this;
        }

        T dot(const Vector3 &rhs) const noexcept {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        T abs_dot(const Vector3 &rhs) const noexcept {
            return std::abs(x * rhs.x + y * rhs.y + z * rhs.z);
        }

        Real angle(const Vector3 &rhs) const noexcept {
            if (dot(rhs) < 0) {
                return pi - 2 * safe_asin(length(*this + rhs) / 2);
            } else {
                return 2 * safe_asin(length(rhs - *this) / 2);
            }
        }

        Vector3 &gram_schmidt(const Vector3 &rhs) noexcept {
            *this -= rhs * dot(rhs);
            return *this;
        }

        Vector3 &cross(const Vector3 &rhs) noexcept {
            T temp_x = difference_of_products(y, rhs.z, z, rhs.y);
            T temp_y = difference_of_products(z, rhs.x, x, rhs.z);
            T temp_z = difference_of_products(x, rhs.y, y, rhs.x);
            x = temp_x;
            y = temp_y;
            z = temp_z;
            return *this;
        }
    };

    template <typename T>
    inline T length_squa(const Vector3<T> &rhs) {
        return rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z;
    }

    template <typename T>
    inline tuple_length_t<T> length(const Vector3<T> &rhs) {
        return std::sqrt(rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z);
    }

    template <typename T>
    inline Vector3<T> normalize(const Vector3<T> &rhs) {
        auto inv_length = tuple_length_t<T> { 1 } / rhs.length();
        return { rhs.x * inv_length, rhs.y * inv_length, rhs.z * inv_length };
    }

    template <typename T>
    inline T dot(const Vector3<T> &lhs, const Vector3<T> &rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    template <typename T>
    inline T abs_dot(const Vector3<T> &lhs, const Vector3<T> &rhs) {
        return std::abs(lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
    }

    template <typename T>
    inline T angle(const Vector3<T> &lhs, const Vector3<T> &rhs) {
        if (dot(lhs, rhs) < 0) {
            return pi - 2 * safe_asin(length(lhs + rhs) / 2);
        } else {
            return 2 * safe_asin(length(rhs - lhs) / 2);
        }
    }

    template <typename T>
    inline Vector3<T> gram_schmidt(const Vector3<T> &lhs, const Vector3<T> &rhs) {
        return lhs - rhs * dot(lhs, rhs);
    }

    template <typename T>
    inline Vector3<T> cross(const Vector3<T> &lhs, const Vector3<T> &rhs) {
        return {
            difference_of_products(lhs.y, rhs.z, lhs.z, rhs.y),
            difference_of_products(lhs.z, rhs.x, lhs.x, rhs.z),
            difference_of_products(lhs.x, rhs.y, lhs.y, rhs.x)
        };
    }

    using Vector2i = Vector2<Int>;
    using Vector2f = Vector2<Real>;
    using Vector3i = Vector3<Int>;
    using Vector3f = Vector3<Real>;

    struct Vector3fi : public  Vector3<Interval> {
        using Vector3<Interval>::x;
        using Vector3<Interval>::y;
        using Vector3<Interval>::z;
        using Vector3<Interval>::has_nan;
        using Vector3<Interval>::operator+;
        using Vector3<Interval>::operator+=;
        using Vector3<Interval>::operator*;
        using Vector3<Interval>::operator*=;

        Vector3fi() noexcept = default;

        Vector3fi(Real x, Real y, Real z) noexcept : Vector3<Interval>(Interval(x), Interval(y), Interval(z)) {}

        Vector3fi(const Interval &x, const Interval &y, const Interval &z) noexcept : Vector3<Interval>(x, y, z) {}

        Vector3fi(const Vector3f &rhs) noexcept : Vector3<Interval>(Interval(rhs.x), Interval(rhs.y), Interval(rhs.z)) {}

        template <typename T>
        explicit Vector3fi(const Point3<T> &rhs) noexcept : Vector3<Interval>(Interval(rhs.x), Interval(rhs.y), Interval(rhs.z)) {}

        Vector3fi(const Vector3<Interval> &rhs) : Vector3<Interval>(rhs) {}

        Vector3fi(const Vector3f &value, const Vector3f &error) noexcept : Vector3<Interval>(
                    Interval::from_value_and_error(value.x, error.x),
                    Interval::from_value_and_error(value.y, error.y),
                    Interval::from_value_and_error(value.z, error.z)
                ) {}

        Vector3f error() const noexcept {
            return { x.width() / 2, y.width() / 2, z.width() / 2 };
        }
        
        bool is_exact() const noexcept {
            return x.width() == 0 && y.width() == 0 && z.width() == 0;
        }
    };
}

template <typename T>
struct fmt::formatter<MatchOfflineRenderer::Vector2<T>> : fmt::formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector2<T> &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector2<T> " << "{ " << v.x << ", " << v.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::Vector2i> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector2i &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector2i " << "{ " << v.x << ", " << v.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::Vector2f> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector2f &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector2f " << "{ " << v.x << ", " << v.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <typename T>
struct fmt::formatter<MatchOfflineRenderer::Vector3<T>> : fmt::formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector3<T> &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector3<T> " << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::Vector3i> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector3i &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector3i " << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::Vector3f> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::Vector3f &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Vector3f " << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};
