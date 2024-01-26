#pragma once

#include <MatchOfflineRenderer/math/point.hpp>

namespace MatchOfflineRenderer {
    // 2D包围盒和3D包围盒
    
    template <typename T>
    struct Bounds2 {
        Point2<T> min, max;

        Bounds2() noexcept {
            T min_num = std::numeric_limits<T>::lowest();
            T max_num = std::numeric_limits<T>::max();
            min = Point2 { min_num, min_num };
            max = Point2 { max_num, max_num };
        }
        
        Bounds2(const Point2<T> &rhs) noexcept : min(rhs), max(rhs) {}
        
        Bounds2(const Point2<T> &lhs, const Point2<T> &rhs) noexcept : min(::MatchOfflineRenderer::min(lhs, rhs)), max(::MatchOfflineRenderer::max(lhs, rhs)) {}

        const Point2<T> &operator[](int i) const noexcept {
            if (i == 0) return min;
            if (i == 1) return max;
        }

        Point2<T> &operator[](int i) noexcept {
            if (i == 0) return min;
            if (i == 1) return max;
        }

        bool is_overlaps(const Bounds2 &rhs) const noexcept {
            return (
                (rhs.max.x >= min.x) && (rhs.min.x <= max.x) &&
                (rhs.max.y >= min.y) && (rhs.min.y <= max.y)
            );
        }

        bool is_contains(const Point2<T> &rhs) const noexcept {
            return (
                (rhs.x >= min.x) && (rhs.x <= max.x) &&
                (rhs.y >= min.y) && (rhs.y <= max.y)
            );
        }

        bool is_contains_exclusive(const Point2<T> &rhs) const noexcept {
            return (
                (rhs.x >= min.x) && (rhs.x < max.x) &&
                (rhs.y >= min.y) && (rhs.y < max.y)
            );
        }

        template <typename R>
        Bounds2 &expand(R delta) noexcept {
            min -= Vector2<T> { delta, delta };
            max += Vector2<T> { delta, delta };
            return *this;
        }

        Vector2<T> diagonal() const noexcept {
            return max - min;
        }

        Real area() const noexcept {
            auto d = diagonal();
            return d.x * d.y;
        }

        int maximum_extent() const noexcept {
            return diagonal().max_dimension();
        }

        Point2<T> lerp(const Point2f &t) const noexcept {
            return Point2 { ::MatchOfflineRenderer::lerp(min.x, max.x, t.x), ::MatchOfflineRenderer::lerp(min.y, max.y, t.y) };
        }
        
        Vector2<T> offset(const Point2<T> &rhs) const noexcept {
            auto v = rhs - min;
            if (max.x > min.x) v.x /= max.x - min.x;
            if (max.y > min.y) v.y /= max.y - min.y;
            return v;
        }

        bool is_empty() const noexcept {
            return min.x >= max.x || min.y >= max.y;
        }

        bool is_degenerate() const noexcept {
            return min.x > max.x || min.y > max.y;
        }

        bool is_inside(const Point2<T> &rhs) const noexcept {
            return (rhs.x >= min.x && rhs.x <= max.x) &&
                (rhs.y >= min.y && rhs.y <= max.y);
        }

        bool is_inside_exclusive(const Point2<T> &rhs) const noexcept {
            return (rhs.x >= min.x && rhs.x < max.x) &&
                (rhs.y >= min.y && rhs.y < max.y);
        }

        inline static Bounds2 generate_union(const Bounds2 &lhs, const Point2<T> rhs) noexcept {
            return {
                ::MatchOfflineRenderer::min(lhs.min, rhs),
                ::MatchOfflineRenderer::max(lhs.max, rhs),
            };
        }

        inline static Bounds2 generate_union(const Bounds2 &lhs, const Bounds2 rhs) noexcept {
            return {
                ::MatchOfflineRenderer::min(lhs.min, rhs.min),
                ::MatchOfflineRenderer::max(lhs.max, rhs.max),
            };
        }

        inline static Bounds2 generate_intersect(const Bounds2 &lhs, const Bounds2 rhs) noexcept {
            return {
                ::MatchOfflineRenderer::max(lhs.min, rhs.min),
                ::MatchOfflineRenderer::min(lhs.max, rhs.max),
            };
        }
    };

    using Bounds2i = Bounds2<Int>;
    using Bounds2f = Bounds2<Real>;

    template <typename T>
    struct Bounds3 {
        Point3<T> min, max;

        struct BoundsSphere {
            Point3<T> center;
            Real radius;
        };

        Bounds3() noexcept {
            T min_num = std::numeric_limits<T>::lowest();
            T max_num = std::numeric_limits<T>::max();
            min = Point3 { min_num, min_num, min_num };
            max = Point3 { max_num, max_num, max_num };
        }
        
        Bounds3(const Point3<T> &rhs) noexcept : min(rhs), max(rhs) {}
        
        Bounds3(const Point3<T> &lhs, const Point3<T> &rhs) noexcept : min(::MatchOfflineRenderer::min(lhs, rhs)), max(::MatchOfflineRenderer::max(lhs, rhs)) {}

        const Point3<T> &operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 2))
            if (i == 0) return min;
            return max;
        }

        Point3<T> &operator[](int i) noexcept {
            MCH_DASSERT((0 <= i) && (i < 2))
            if (i == 0) return min;
            return max;
        }

        Point3<T> corner(int i) const noexcept {
            return Point3 {
                (*this)[i & 1].x,
                (*this)[(i & 2) ? 1 : 0].y,
                (*this)[(i & 4) ? 1 : 0 ].z
            };
        }

        bool is_overlaps(const Bounds3 &rhs) const noexcept {
            return (
                (rhs.max.x >= min.x) && (rhs.min.x <= max.x) &&
                (rhs.max.y >= min.y) && (rhs.min.y <= max.y) &&
                (rhs.max.z >= min.z) && (rhs.min.z <= max.z)
            );
        }

        bool is_contains(const Point3<T> &rhs) const noexcept {
            return (
                (rhs.x >= min.x) && (rhs.x <= max.x) &&
                (rhs.y >= min.y) && (rhs.y <= max.y) &&
                (rhs.z >= min.z) && (rhs.z <= max.z)
            );
        }

        bool is_contains_exclusive(const Point3<T> &rhs) const noexcept {
            return (
                (rhs.x >= min.x) && (rhs.x < max.x) &&
                (rhs.y >= min.y) && (rhs.y < max.y) &&
                (rhs.z >= min.z) && (rhs.z < max.z)
            );
        }

        template <typename R>
        Bounds3 &expand(R delta) noexcept {
            min -= Vector3<T> { delta, delta, delta };
            max += Vector3<T> { delta, delta, delta };
            return *this;
        }

        Vector3<T> diagonal() const noexcept {
            return max - min;
        }

        Real area() const noexcept {
            auto d = diagonal();
            return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
        }

        Real volume() const noexcept {
            auto d = diagonal();
            return d.x * d.y * d.z;
        }

        int maximum_extent() const noexcept {
            return diagonal().max_dimension();
        }

        Point3<T> lerp(const Point3f &t) const noexcept {
            return Point3 { lerp(min.x, max.x, t.x), lerp(min.y, max.y, t.y), lerp(min.z, max.z, t.z) };
        }
        
        Vector3<T> offset(const Point3<T> &rhs) const noexcept {
            auto v = rhs - min;
            if (max.x > min.x) v.x /= max.x - min.x;
            if (max.y > min.y) v.y /= max.y - min.y;
            if (max.z > min.z) v.z /= max.z - min.z;
            return v;
        }

        bool is_empty() const noexcept {
            return min.x >= max.x || min.y >= max.y || min.z >= max.z;
        }

        bool is_degenerate() const noexcept {
            return min.x > max.x || min.y > max.y || min.z > max.z;
        }

        bool is_inside(const Point3<T> &rhs) const noexcept {
            return (rhs.x >= min.x && rhs.x <= max.x) &&
                (rhs.y >= min.y && rhs.y <= max.y) &&
                (rhs.z >= min.z && rhs.z <= max.z);
        }

        bool is_inside_exclusive(const Point3<T> &rhs) const noexcept {
            return (rhs.x >= min.x && rhs.x < max.x) &&
                (rhs.y >= min.y && rhs.y < max.y) &&
                (rhs.z >= min.z && rhs.z < max.z);
        }

        BoundsSphere generate_bounds_sphere() const noexcept {
            auto center = (min + max) / 2;
            return { center, is_inside(center) ? distance(center, max) : 0 };
        }

        inline static Bounds3 generate_union(const Bounds3 &lhs, const Point3<T> rhs) noexcept {
            return {
                ::MatchOfflineRenderer::min(lhs.min, rhs),
                ::MatchOfflineRenderer::max(lhs.max, rhs),
            };
        }

        inline static Bounds3 generate_union(const Bounds3 &lhs, const Bounds3 rhs) noexcept {
            return {
                ::MatchOfflineRenderer::min(lhs.min, rhs.min),
                ::MatchOfflineRenderer::max(lhs.max, rhs.max),
            };
        }

        inline static Bounds3 generate_intersect(const Bounds3 &lhs, const Bounds3 rhs) noexcept {
            return {
                ::MatchOfflineRenderer::max(lhs.min, rhs.min),
                ::MatchOfflineRenderer::min(lhs.max, rhs.max),
            };
        }
    };

    using Bounds3i = Bounds3<Int>;
    using Bounds3f = Bounds3<Real>;

    class Bounds2iIterator : public std::forward_iterator_tag {
    public:
        Bounds2iIterator(const Bounds2i &bounds, const Point2i &current)
            : current(current), bounds(&bounds) {}

        Bounds2iIterator operator++() {
            advance();
            return *this;
        }

        Bounds2iIterator operator++(int) {
            Bounds2iIterator old = *this;
            advance();
            return old;
        }

        bool operator==(const Bounds2iIterator &rhs) const {
            return current == rhs.current && bounds == rhs.bounds;
        }

        bool operator!=(const Bounds2iIterator &rhs) const {
            return current != rhs.current || bounds != rhs.bounds;
        }

        Point2i operator*() const { return Point2i { current }; }

    private:
        void advance() {
            ++current.x;
            if (current.x == bounds->max.x) {
                current.x = bounds->min.x;
                ++current.y;
            }
        }
        
        Point2i current;
        const Bounds2i *bounds;
    };

    inline Bounds2iIterator begin(const Bounds2i &bounds) {
        return Bounds2iIterator(bounds, bounds.min);
    }

    inline Bounds2iIterator end(const Bounds2i &bounds) {
        if (bounds.min.x >= bounds.max.x || bounds.min.y >= bounds.max.y) {
            return Bounds2iIterator(bounds, bounds.min);
        }
        return Bounds2iIterator(bounds, { bounds.min.x, bounds.max.y });
    }
}
