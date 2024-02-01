#pragma once

#include <MatchOfflineRenderer/math/vector.hpp>

namespace MatchOfflineRenderer::math {
    // 2D点和3D点

    template <typename T>
    struct Point2 : public Tuple2<Point2, T> {
        using Tuple2<Point2, T>::x;
        using Tuple2<Point2, T>::y;
        using Tuple2<Point2, T>::has_nan;
        using Tuple2<Point2, T>::operator+;
        using Tuple2<Point2, T>::operator+=;
        using Tuple2<Point2, T>::operator*;
        using Tuple2<Point2, T>::operator*=;

        Point2() noexcept = default;
       
        Point2(T x, T y) noexcept : Tuple2<Point2, T>(x, y) {}

        template <typename R>
        explicit Point2(const Point2<R> &rhs) noexcept : Tuple2<Point2, T>(static_cast<T>(rhs.x), static_cast<T>(rhs.y)) {}
        
        template <typename R>
        explicit Point2(const Vector2<R> &rhs) noexcept : Tuple2<Point2, T>(static_cast<T>(rhs.x), static_cast<T>(rhs.y)) {}

        Point2<T> operator-() const noexcept { return { -x, -y }; };

        template <typename R>
        Point2<decltype(T {} + R {})> operator+(const Vector2<R> &rhs) const noexcept {
            return { x + rhs.x, y + rhs.y };
        }
        
        template <typename R>
        Point2<T> &operator+=(const Vector2<R> &rhs) const noexcept {
            x += rhs.x;
            y += rhs.y;
            return *this;
        }

        template <typename R>
        Point2<decltype(T {} - R {})> operator-(const Vector2<R> &rhs) const noexcept {
            return { x - rhs.x, y - rhs.y };
        }
        
        template <typename R>
        Point2<T> &operator-=(const Vector2<R> &rhs) const noexcept {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }
        
        template <typename R>
        Vector2<decltype(T {} - R {})> operator-(const Point2<R> &rhs) const noexcept {
            return { x - rhs.x, y - rhs.y };
        }

        T distance_squa(const Point2<T> &rhs) const noexcept {
            return (*this - rhs).length_squa();
        }

        tuple_length_t<T> distance(const Point2<T> &rhs) const noexcept {
            return (*this - rhs).length();
        }
    };

    template <typename T>
    inline T distance_squa(const Point2<T> &lhs, const Point2<T> &rhs) noexcept {
        return (lhs - rhs).length_squa();
    }

    template <typename T>
    inline tuple_length_t<T> distance(const Point2<T> &lhs, const Point2<T> &rhs) noexcept {
        return (lhs - rhs).length();
    }

    template <typename T>
    struct Point3 : public Tuple3<Point3, T> {
        using Tuple3<Point3, T>::x;
        using Tuple3<Point3, T>::y;
        using Tuple3<Point3, T>::z;
        using Tuple3<Point3, T>::has_nan;
        using Tuple3<Point3, T>::operator+;
        using Tuple3<Point3, T>::operator+=;
        using Tuple3<Point3, T>::operator*;
        using Tuple3<Point3, T>::operator*=;

        Point3() noexcept = default;
       
        Point3(T x, T y, T z) noexcept : Tuple3<Point3, T>(x, y, z) {}

        template <typename R>
        explicit Point3(const Point3<R> &rhs) noexcept : Tuple3<Point3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}
        
        template <typename R>
        explicit Point3(const Vector3<R> &rhs) noexcept : Tuple3<Point3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

        Point3<T> operator-() const noexcept { return { -x, -y, -z }; };

        template <typename R>
        Point3<decltype(T {} + R {})> operator+(const Vector3<R> &rhs) const noexcept {
            return { x + rhs.x, y + rhs.y, z + rhs.z };
        }
        
        template <typename R>
        Point3<T> &operator+=(const Vector3<R> &rhs) const noexcept {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        template <typename R>
        Point3<decltype(T {} - R {})> operator-(const Vector3<R> &rhs) const noexcept {
            return { x - rhs.x, y - rhs.y, z - rhs.z };
        }
        
        template <typename R>
        Point3<T> &operator-=(const Vector3<R> &rhs) const noexcept {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }
        
        template <typename R>
        Vector3<decltype(T {} - R {})> operator-(const Point3<R> &rhs) const noexcept {
            return { x - rhs.x, y - rhs.y, z - rhs.z };
        }

        T distance_squa(const Point3<T> &rhs) const noexcept {
            return (*this - rhs).length_squa();
        }

        tuple_length_t<T> distance(const Point3<T> &rhs) const noexcept {
            return (*this - rhs).length();
        }
    };

    template <typename T>
    inline T distance_squa(const Point3<T> &lhs, const Point3<T> &rhs) noexcept {
        return (lhs - rhs).length_squa();
    }

    template <typename T>
    inline tuple_length_t<T> distance(const Point3<T> &lhs, const Point3<T> &rhs) noexcept {
        return (lhs - rhs).length();
    }

    template <typename T>
    template <typename R>
    Vector2<T>::Vector2(const Point2<R> &rhs) noexcept : Tuple3<Vector2, T>(T { rhs.x }, T { rhs.y }) {}

    template <typename T>
    template <typename R>
    Vector3<T>::Vector3(const Point3<R> &rhs) noexcept : Tuple3<Vector3, T>(T { rhs.x }, T { rhs.y }, T { rhs.z }) {}

    using Point2i = Point2<Int>;
    using Point2f = Point2<Real>;
    using Point3i = Point3<Int>;
    using Point3f = Point3<Real>;


    struct Point3fi : public  Point3<Interval> {
        using Point3<Interval>::x;
        using Point3<Interval>::y;
        using Point3<Interval>::z;
        using Point3<Interval>::has_nan;
        using Point3<Interval>::operator+;
        using Point3<Interval>::operator+=;
        using Point3<Interval>::operator*;
        using Point3<Interval>::operator*=;

        Point3fi() noexcept = default;

        Point3fi(Real x, Real y, Real z) noexcept : Point3<Interval>(Interval(x), Interval(y), Interval(z)) {}

        Point3fi(const Interval &x, const Interval &y, const Interval &z) noexcept : Point3<Interval>(x, y, z) {}

        Point3fi(const Vector3f &rhs) noexcept : Point3<Interval>(Interval(rhs.x), Interval(rhs.y), Interval(rhs.z)) {}

        template <typename T>
        explicit Point3fi(const Point3<T> &rhs) noexcept : Point3<Interval>(Interval(rhs.x), Interval(rhs.y), Interval(rhs.z)) {}

        Point3fi(const Point3<Interval> &rhs) : Point3<Interval>(rhs) {}

        Point3fi(const Vector3f &value, const Vector3f &error) noexcept : Point3<Interval>(
                    Interval::from_value_and_error(value.x, error.x),
                    Interval::from_value_and_error(value.y, error.y),
                    Interval::from_value_and_error(value.z, error.z)
                ) {}

        Point3f error() const noexcept {
            return { x.width() / 2, y.width() / 2, z.width() / 2 };
        }
        
        bool is_exact() const noexcept {
            return x.width() == 0 && y.width() == 0 && z.width() == 0;
        }
    };
}

template <typename T>
struct fmt::formatter<MatchOfflineRenderer::math::Point2<T>> : fmt::formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point2<T> &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point2<T> " << "{ " << p.x << ", " << p.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::math::Point2i> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point2i &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point2i " << "{ " << p.x << ", " << p.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::math::Point2f> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point2f &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point2f " << "{ " << p.x << ", " << p.y << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <typename T>
struct fmt::formatter<MatchOfflineRenderer::math::Point3<T>> : fmt::formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point3<T> &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point3<T> " << "{ " << p.x << ", " << p.y << ", " << p.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::math::Point3i> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point3i &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point3i " << "{ " << p.x << ", " << p.y << ", " << p.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};

template <>
struct fmt::formatter<MatchOfflineRenderer::math::Point3f> : formatter<std::string_view> {
    auto format(const MatchOfflineRenderer::math::Point3f &p, format_context& ctx) const noexcept {
        std::stringstream ss;
        ss << "Point3f " << "{ " << p.x << ", " << p.y << ", " << p.z << " }";
        return formatter<string_view>::format(ss.str(), ctx);
    }
};
