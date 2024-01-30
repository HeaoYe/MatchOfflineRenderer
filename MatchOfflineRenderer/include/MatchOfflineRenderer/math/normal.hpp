#pragma once

#include <MatchOfflineRenderer/math/vector.hpp>

namespace MatchOfflineRenderer::math {
    // 3D法线
    
    template <typename T>
    struct Normal3 : public Tuple3<Normal3, T> {
        using Tuple3<Normal3, T>::x;
        using Tuple3<Normal3, T>::y;
        using Tuple3<Normal3, T>::z;

        Normal3() noexcept = default;

        Normal3(T x, T y, T z) noexcept : Tuple3<Normal3, T>(x, y, z) {}

        template <typename R>
        explicit Normal3(const Normal3<R> &rhs) noexcept : Tuple3<Normal3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

        template <typename R>
        explicit Normal3(const Vector3<R> &rhs) noexcept : Tuple3<Normal3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

        T length_squa() const noexcept {
            return x * x + y * y + z * z;
        }

        tuple_length_t<T> length() const noexcept {
            return std::sqrt<tuple_length_t<T>>(x * x + y * y + z * z);
        }

        Normal3 &normalize() noexcept {
            auto inv_length = tuple_length_t<T> { 1 } / length();
            x *= inv_length;
            y *= inv_length;
            z *= inv_length;
            return *this;
        }

        T dot(const Normal3 &rhs) const noexcept {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        T dot(const Vector3<T> &rhs) const noexcept {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        T abs_dot(const Normal3 &rhs) const noexcept {
            return std::abs(x * rhs.x + y * rhs.y + z * rhs.z);
        }

        T abs_dot(const Vector3<T> &rhs) const noexcept {
            return std::abs(x * rhs.x + y * rhs.y + z * rhs.z);
        }

        Real angle(const Normal3 &rhs) const noexcept {
            if (dot(rhs) < 0) {
                return pi - 2 * safe_asin(length(*this + rhs) / 2);
            } else {
                return 2 * safe_asin(length(rhs - *this) / 2);
            }
        }

        Normal3 gram_schmidt(const Normal3 &rhs) {
            return *this - rhs * dot(rhs);
        }

        Normal3 cross(const Normal3 &rhs) {
            return {
                difference_of_products(y, rhs.z, z, rhs.y),
                difference_of_products(z, rhs.x, x, rhs.z),
                difference_of_products(x, rhs.y, y, rhs.x)
            };
        }
    };

    template <typename T>
    inline T length_squa(const Normal3<T> &rhs) {
        return rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z;
    }

    template <typename T>
    inline tuple_length_t<T> length(const Normal3<T> &rhs) {
        return std::sqrt<tuple_length_t<T>>(rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z);
    }

    template <typename T>
    inline Normal3<T> normalize(const Normal3<T> &rhs) {
        auto inv_length = tuple_length_t<T> { 1 } / rhs.length();
        return { rhs.x * inv_length, rhs.y * inv_length, rhs.z * inv_length };
    }

    template <typename T>
    inline T dot(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    template <typename T>
    inline T dot(const Normal3<T> &lhs, const Vector3<T> &rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    template <typename T>
    inline T abs_dot(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        return std::abs(lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
    }

    template <typename T>
    inline T abs_dot(const Normal3<T> &lhs, const Vector3<T> &rhs) {
        return std::abs(lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
    }

    template <typename T>
    inline T angle(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        if (dot(lhs, rhs) < 0) {
            return pi - 2 * safe_asin(length(lhs + rhs) / 2);
        } else {
            return 2 * safe_asin(length(rhs - lhs) / 2);
        }
    }

    template <typename T>
    inline Normal3<T> gram_schmidt(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        return lhs - rhs * dot(lhs, rhs);
    }

    template <typename T>
    inline Normal3<T> cross(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        return {
            difference_of_products(lhs.y, rhs.z, lhs.z, rhs.y),
            difference_of_products(lhs.z, rhs.x, lhs.x, rhs.z),
            difference_of_products(lhs.x, rhs.y, lhs.y, rhs.x)
        };
    }

    template <typename T>
    inline Normal3<T> face_forward(const Normal3<T> &lhs, const Normal3<T> &rhs) {
        return (dot(lhs, rhs) < 0.0f) ? -lhs : lhs;
    }

    template <typename T>
    template <typename R>
    Vector3<T>::Vector3(const Normal3<R> &rhs) noexcept : Tuple3<Vector3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

    using Normal3f = Normal3<Real>;
}

template <typename T>
struct fmt::formatter<MatchOfflineRenderer::math::Normal3<T>> : fmt::formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Normal3<T> &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Normal3<T> " << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::math::Normal3f> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Normal3f &v, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Normal3f " << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};
