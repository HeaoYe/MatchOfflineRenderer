#pragma once

#include <MatchOfflineRenderer/math/types.hpp>

#include <cmath>

namespace MatchOfflineRenderer::math {
    // 2元数和3元数,是 向量 点 法线 的基类
    template <template <typename> class Child, typename T>
    struct Tuple2 {
        static constexpr int nDimensions = 2;
        T x {}, y {};

        Tuple2() noexcept = default;
        Tuple2(T x, T y) noexcept : x(x), y(y) {
            MCH_DASSERT(!has_nan())
        }

        bool has_nan() const noexcept {
            return isnan(x) || isnan(y);
        }

        T operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 2))
            if (i == 0) return x;
            return y;
        }

        T &operator[](int i) noexcept {
            MCH_DASSERT((0 <= i) && (i < 2))
            if (i == 0) return x;
            return y;
        }

        bool operator==(const Child<T> &rhs) const noexcept {
            return (x == rhs.x) && (y == rhs.y);
        }

        bool operator!=(const Child<T> &rhs) const noexcept {
            return (x != rhs.x) || (y != rhs.y);
        }

        template <typename R>
        Child<decltype(T {} + R {})> operator+(const Child<R> &rhs) const noexcept {
            MCH_DASSERT(!rhs.has_nan())
            return { x + rhs.x, y + rhs.y };
        }

        template <typename R>
        Child<T> &operator+=(const Child<R> &rhs) noexcept {
            MCH_DASSERT(!rhs.has_nan())
            x += rhs.x;
            y += rhs.y;
            return static_cast<Child<T> &>(*this);
        }

        template <typename R>
        Child<decltype(T {} - R {})> operator-(const Child<R> &rhs) const noexcept {
            MCH_DASSERT(!rhs.has_nan())
            return { x - rhs.x, y - rhs.y };
        }

        template <typename R>
        Child<T> &operator-=(const Child<R> &rhs) noexcept {
            MCH_DASSERT(!rhs.has_nan())
            x -= rhs.x;
            y -= rhs.y;
            return static_cast<Child<T> &>(*this);
        }

        Child<T> operator-() const noexcept {
            return Child<T> { -x, -y };
        }

        template <typename R>
        Child<decltype(T {} * R {})> operator*(R rhs) const noexcept {
            MCH_DASSERT(!isnan(rhs))
            return { x * rhs, y * rhs };
        }

        template <typename R>
        Child<T> &operator*=(R rhs) noexcept {
            MCH_DASSERT(!isnan(rhs))
            x *= rhs;
            y *= rhs;
            return static_cast<Child<T> &>(*this);
        }

        template <typename R>
        Child<decltype(T {} / R {})> operator/(R rhs) const noexcept {
            rhs = static_cast<R>(1) / rhs;
            MCH_DASSERT(!isnan(rhs))
            return { x * rhs, y * rhs };
        }

        template <typename R>
        Child<T> &operator/=(R rhs) noexcept {
            rhs = static_cast<R>(1) / rhs;
            MCH_DASSERT(!isnan(rhs))
            x *= rhs;
            y *= rhs;
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &abs() noexcept {
            x = std::abs(x);
            y = std::abs(y);
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &ceil() noexcept {
            x = std::ceil(x);
            y = std::ceil(y);
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &floor() noexcept {
            x = std::floor(x);
            y = std::floor(y);
            return static_cast<Child<T> &>(*this);
        }

        T min() const noexcept {
            return std::min(x, y);
        }

        T max() const noexcept {
            return std::max(x, y);
        }

        T min_dimension() const noexcept {
            return (x < y) ? 0 : 1;
        }

        T max_dimension() const noexcept {
            return (x > y) ? 0 : 1;
        }

        T h_prod() const noexcept {
            return x * y;
        }
    };

    template <template <typename> class Child, typename T, typename R>
    inline Child<decltype(T {} * R {})> operator*(R lhs, const Tuple2<Child, T> &rhs) noexcept {
        MCH_DASSERT(!isnan(lhs))
        MCH_DASSERT(!rhs.has_nan())
        return rhs * lhs;
    }

    template <template <typename> class Child, typename T>
    inline Child<T> abs(const Tuple2<Child, T> &rhs) noexcept {
        return { std::abs(rhs.x), std::abs(rhs.y) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> ceil(const Tuple2<Child, T> &rhs) noexcept {
        return { std::ceil(rhs.x), std::ceil(rhs.y) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> floor(const Tuple2<Child, T> &rhs) noexcept {
        return { std::floor(rhs.x), std::floor(rhs.y) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> lerp(const Tuple2<Child, T> &lhs, const Tuple2<Child, T> &rhs, Real t) noexcept {
        return (1 - t) * lhs + t * rhs;
    }

    template <template <typename> class Child, typename T>
    inline Child<T> fma(Real a, const Tuple2<Child, T> &b, const Tuple2<Child, T> &c) noexcept {
        return { a * b.x + c.x, a * b.y + c.y };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> min(const Tuple2<Child, T> &lhs, const Tuple2<Child, T> &rhs) noexcept {
        return { std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> max(const Tuple2<Child, T> &lhs, const Tuple2<Child, T> &rhs) noexcept {
        return { std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y) };
    }

    template <template <typename> class Child, typename T>
    inline T min_value(const Tuple2<Child, T> &rhs) noexcept {
        return rhs.min();
    }

    template <template <typename> class Child, typename T>
    inline T max_value(const Tuple2<Child, T> &rhs) noexcept {
        return rhs.max();
    }

    template <template <typename> class Child, typename T>
    inline T min_dimension(const Tuple2<Child, T> &rhs) noexcept {
        return rhs.min_dimension();
    }

    template <template <typename> class Child, typename T>
    inline T max_dimension(const Tuple2<Child, T> &rhs) noexcept {
        return rhs.max_dimension();
    }

    template <template <typename> class Child, typename T>
    inline Child<T> permute(const Tuple2<Child, T> &rhs, int x, int y) noexcept {
        return { rhs[x], rhs[y] };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> h_prod(const Tuple2<Child, T> &rhs) noexcept {
        return rhs.h_prod();
    }

    template <template <typename> class Child, typename T>
    struct Tuple3 {
        static constexpr int nDimensions = 3;
        T x {}, y {}, z {};

        Tuple3() noexcept = default;
        Tuple3(T x, T y, T z) noexcept : x(x), y(y), z(z) {
            MCH_DASSERT(!has_nan())
        }

        bool has_nan() const noexcept {
            return isnan(x) || isnan(y) || isnan(z);
        }

        T operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        T &operator[](int i) noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        bool operator==(const Child<T> &rhs) const noexcept {
            return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);
        }

        bool operator!=(const Child<T> &rhs) const noexcept {
            return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
        }

        template <typename R>
        Child<decltype(T {} + R {})> operator+(const Child<R> &rhs) const noexcept {
            MCH_DASSERT(!rhs.has_nan())
            return { x + rhs.x, y + rhs.y, z + rhs.z };
        }

        template <typename R>
        Child<T> &operator+=(const Child<R> &rhs) noexcept {
            MCH_DASSERT(!rhs.has_nan())
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return static_cast<Child<T> &>(*this);
        }

        template <typename R>
        Child<decltype(T {} - R {})> operator-(const Child<R> &rhs) const noexcept {
            MCH_DASSERT(!rhs.has_nan())
            return { x - rhs.x, y - rhs.y, z - rhs.z };
        }

        template <typename R>
        Child<T> &operator-=(const Child<R> &rhs) noexcept {
            MCH_DASSERT(!rhs.has_nan())
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return static_cast<Child<T> &>(*this);
        }

        Child<T> operator-() const noexcept {
            return Child<T> { -x, -y, -z };
        }

        template <typename R>
        Child<decltype(T {} * R {})> operator*(R rhs) const noexcept {
            MCH_DASSERT(!isnan(rhs))
            return { x * rhs, y * rhs, z * rhs };
        }

        template <typename R>
        Child<T> &operator*=(R rhs) noexcept {
            MCH_DASSERT(!isnan(rhs))
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return static_cast<Child<T> &>(*this);
        }

        template <typename R>
        Child<decltype(T {} / R {})> operator/(R rhs) const noexcept {
            rhs = static_cast<R>(1) / rhs;
            MCH_DASSERT(!isnan(rhs))
            return { x * rhs, y * rhs, z * rhs };
        }

        template <typename R>
        Child<T> &operator/=(R rhs) noexcept {
            rhs = static_cast<R>(1) / rhs;
            MCH_DASSERT(!isnan(rhs))
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &abs() noexcept {
            x = std::abs(x);
            y = std::abs(y);
            z = std::abs(z);
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &ceil() noexcept {
            x = std::ceil(x);
            y = std::ceil(y);
            z = std::ceil(z);
            return static_cast<Child<T> &>(*this);
        }

        Child<T> &floor() noexcept {
            x = std::floor(x);
            y = std::floor(y);
            z = std::floor(z);
            return static_cast<Child<T> &>(*this);
        }

        T min() const noexcept {
            return std::min(std::min(x, y), z);
        }

        T max() const noexcept {
            return std::max(std::max(x, y), z);
        }

        Int min_dimension() const noexcept {
            return (x < y) ? ((x < z) ? 0 : 2) : ((y < z) ? 1 : 2);
        }

        Int max_dimension() const noexcept {
            return (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
        }

        T h_prod() const noexcept {
            return x * y * z;
        }
    };

    template <template <typename> class Child, typename T, typename R>
    inline Child<decltype(T {} * R {})> operator*(R lhs, const Tuple3<Child, T> &rhs) noexcept {
        MCH_DASSERT(!isnan(lhs))
        MCH_DASSERT(!rhs.has_nan())
        return rhs * lhs;
    }

    template <template <typename> class Child, typename T>
    inline Child<T> abs(const Tuple3<Child, T> &rhs) noexcept {
        return { std::abs(rhs.x), std::abs(rhs.y), std::abs(rhs.z) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> ceil(const Tuple3<Child, T> &rhs) noexcept {
        return { std::ceil(rhs.x), std::ceil(rhs.y), std::ceil(rhs.z) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> floor(const Tuple3<Child, T> &rhs) noexcept {
        return { std::floor(rhs.x), std::floor(rhs.y), std::floor(rhs.z) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> lerp(const Tuple3<Child, T> &lhs, const Tuple3<Child, T> &rhs, Real t) noexcept {
        return lhs * (1 - t) + rhs * t;
    }

    template <template <typename> class Child, typename T>
    inline Child<T> fma(Real a, const Tuple3<Child, T> &b, const Tuple3<Child, T> &c) noexcept {
        return { a * b.x + c.x, a * b.y + c.y, a * b.z + c.z };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> min(const Tuple3<Child, T> &lhs, const Tuple3<Child, T> &rhs) noexcept {
        return { std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y), std::min(lhs.z, rhs.z) };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> max(const Tuple3<Child, T> &lhs, const Tuple3<Child, T> &rhs) noexcept {
        return { std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y), std::max(lhs.z, rhs.z) };
    }

    template <template <typename> class Child, typename T>
    inline T min_value(const Tuple3<Child, T> &rhs) noexcept {
        return rhs.min();
    }

    template <template <typename> class Child, typename T>
    inline T max_value(const Tuple3<Child, T> &rhs) noexcept {
        return rhs.max();
    }

    template <template <typename> class Child, typename T>
    inline T min_dimension(const Tuple3<Child, T> &rhs) noexcept {
        return rhs.min_dimension();
    }

    template <template <typename> class Child, typename T>
    inline T max_dimension(const Tuple3<Child, T> &rhs) noexcept {
        return rhs.max_dimension();
    }

    template <template <typename> class Child, typename T>
    inline Child<T> permute(const Tuple3<Child, T> &rhs, int x, int y, int z) noexcept {
        return { rhs[x], rhs[y], rhs[z] };
    }

    template <template <typename> class Child, typename T>
    inline Child<T> h_prod(const Tuple3<Child, T> &rhs) noexcept {
        return rhs.h_prod();
    }

    template <typename T>
    struct tuple_length {
        using type = Real;
    };

    template <>
    struct tuple_length<double> {
        using type = double;
    };

    template <typename T>
    using tuple_length_t = typename tuple_length<T>::type;
}
