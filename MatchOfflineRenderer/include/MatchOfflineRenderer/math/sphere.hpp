#pragma once

#include <MatchOfflineRenderer/math/vector.hpp>
#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/bounds.hpp>
#include <MatchOfflineRenderer/math/transform.hpp>

namespace MatchOfflineRenderer {
    // 球面参数化
    // 实现了三种参数化
    // 1.球坐标
    // 2.八面体映射, 用两个float存储球面上的某一点,映射前后面积不成比例
    // 3.等积映射(好像叫这个),用两个float存储球面上的某一点,映射前后面积成比例

    inline Real spherical_triangle_area(const Vector3f &a, const Vector3f &b, const Vector3f &c) noexcept {
        return std::abs(2 * std::atan2(dot(a, cross(b, c)), 1 + dot(a, b) + dot(a, c) + dot(b, c)));
    }
    
    inline Real spherical_quad_area(const Vector3f &a, const Vector3f &b, const Vector3f &c, const Vector3f &d) noexcept {
        Vector3f axb = cross(a, b), bxc = cross(b, c);
        Vector3f cxd = cross(c, d), dxa = cross(d, a);
        if (axb.length_squa() == 0 || bxc.length_squa() == 0 || cxd.length_squa() == 0 || dxa.length_squa() == 0)
            return 0;
        axb.normalize();
        bxc.normalize();
        cxd.normalize();
        dxa.normalize();

        Real alpha = angle(dxa, -axb);
        Real beta = angle(axb, -bxc);
        Real gamma = angle(bxc, -cxd);
        Real delta = angle(cxd, -dxa);

        return std::abs(alpha + beta + gamma + delta - 2 * pi);
    }

    inline Vector3f spherical_direction(Real sin_theta, Real cos_theta, Real phi) noexcept {
        return {
            std::clamp<Real>(sin_theta, -1, 1) * std::cos(phi),
            std::clamp<Real>(sin_theta, -1, 1) * std::sin(phi),
            std::clamp<Real>(cos_theta, -1, 1)
        };
    }

    inline Real spherical_theta(const Vector3f &rhs) noexcept {
        return safe_acos(rhs.z);
    }

    inline Real spherical_phi(const Vector3f &rhs) noexcept {
        Real p = std::atan2(rhs.y, rhs.x);
        return (p < 0) ?  (p + 2 * pi) : p;
    }

    inline Real spherical_cos_theta(const Vector3f &rhs) noexcept {
        return rhs.z;
    };

    inline Real spherical_abs_cos_theta(const Vector3f &rhs) noexcept {
        return std::abs(rhs.z);
    }
    
    inline Real spherical_cos2_theta(const Vector3f &rhs) noexcept {
        return rhs.z * rhs.z;
    }

    inline Real spherical_sin2_theta(const Vector3f &rhs) noexcept {
        return std::max<Real>(0, 1 - spherical_cos2_theta(rhs));
    }

    inline Real spherical_sin_theta(const Vector3f &rhs) noexcept {
        return std::sqrt(spherical_sin2_theta(rhs));
    }

    inline Real spherical_tan_theta(const Vector3f &rhs) noexcept {
        return spherical_sin_theta(rhs) / spherical_cos_theta(rhs);
    }

    inline Real spherical_tan2_theta(const Vector3f &rhs) noexcept {
        return spherical_sin2_theta(rhs) / spherical_cos2_theta(rhs);
    }

    inline Real spherical_cos_phi(const Vector3f &rhs) noexcept {
        Real sin_theta = spherical_sin_theta(rhs);
        return (sin_theta == 0) ? 1 : std::clamp<Real>(rhs.x / sin_theta, -1, 1);
    }

    inline Real spherical_sin_phi(const Vector3f &rhs) noexcept {
        Real sin_theta = spherical_sin_theta(rhs);
        return (sin_theta == 0) ? 0 : std::clamp<Real>(rhs.y / sin_theta, -1, 1);
    }

    inline Real spherical_cos_delta_phi(const Vector3f &lhs, const Vector3f &rhs) noexcept {
        Real lhs_xy = lhs.x * lhs.x + lhs.y * lhs.y;
        Real rhs_xy = rhs.x * rhs.x + rhs.y * rhs.y;
        if (lhs_xy == 0 || rhs_xy == 0) {
            return 1;
        }
        return std::clamp<Real>((lhs.x * rhs.x + lhs.y * rhs.y) / std::sqrt(lhs_xy * rhs_xy), -1, 1);
    }

    struct OctahedralVector {
        uint16_t x {}, y {};

        OctahedralVector(const Vector3f &rhs) noexcept {
            auto v = rhs / (std::abs(rhs.x) + std::abs(rhs.y) + std::abs(rhs.z));
            if (v.z >= 0) {
                x = real_encode_uint16(v.x);
                y = real_encode_uint16(v.y);
            } else {
                x = real_encode_uint16((1 - std::abs(v.y)) * std::copysign(1, v.x));
                x = real_encode_uint16((1 - std::abs(v.x)) * std::copysign(1, v.y));
            }
        }

        explicit operator Vector3f() const noexcept {
            Vector3f result {
                -1 + 2 * (static_cast<Real>(x) / Real { 65535 }),
                -1 + 2 * (static_cast<Real>(y) / Real { 65535 }),
                0
            };
            result.z = 1 - (std::abs(result.x) + std::abs(result.y));
            if (result.z < 0) {
                Real temp_x = result.x;
                result.x = (1 - std::abs(result.y)) * std::copysign(1, result.x);
                result.y = (1 - std::abs(temp_x)) * std::copysign(1, result.y);
            }
            return result.normalize();
        }
    };

    inline Vector3f equal_area_square_to_sphere(const Point2f &rhs) noexcept {
        Real u = 2 * rhs.x - 1;
        Real v = 2 * rhs.y - 1;
        Real up = std::abs(u);
        Real vp = std::abs(v);

        Real signed_distance = 1 - (up * vp);
        Real d = std::abs(signed_distance);
        Real r = 1 - d;
        Real r2 = r * r;

        Real z = std::copysign(1 - r2, signed_distance);

        Real phi = (r == 0) ? 1 : (((vp - up) / r) + 1) * pi_div4;
        Real cos_phi = std::copysign(std::cos(phi), u);
        Real sin_phi = std::copysign(std::sin(phi), v);

        return { cos_phi * r * safe_sqrt(2 - r2), sin_phi * r * safe_sqrt(2 - r2), z };
    }

    struct DirectionCone {
        Vector3f w {};
        Real cos_theta = std::numeric_limits<Real>::infinity();

        DirectionCone() noexcept = default;

        DirectionCone(const Vector3f &rhs, Real cos_theta) noexcept : w(normalize(rhs)), cos_theta(cos_theta) {}

        explicit DirectionCone(const Vector3f &rhs) noexcept : DirectionCone(rhs, 1) {}

        inline static DirectionCone generate_entire_sphere() {
            return { Vector3f { 0, 0, 1 }, -1 };
        }

        inline static DirectionCone generate_bounds_subtended_directions(const Bounds3f &bounds, const Point3f &point) {
            auto bounds_sphere = bounds.generate_bounds_sphere();
            bounds_sphere.radius *= bounds_sphere.radius;
            auto v = point - bounds_sphere.center;
            auto vlq = v.length_squa();
            if (vlq < bounds_sphere.radius) {
                return generate_entire_sphere();
            }
            return { normalize(v), safe_sqrt(1 - bounds_sphere.radius / vlq) };
        }

        inline static DirectionCone generate_union(const DirectionCone &lhs, const DirectionCone &rhs) {
            if (lhs.is_empty()) return rhs;
            if (rhs.is_empty()) return lhs;

            Real theta_lhs = safe_acos(lhs.cos_theta);
            Real theta_rhs = safe_acos(rhs.cos_theta);
            Real theta_delta = angle(lhs.w, rhs.w);
            if (std::min(theta_delta + theta_rhs, pi) <= theta_lhs) return lhs;
            if (std::min(theta_delta + theta_lhs, pi) <= theta_rhs) return rhs;

            Real theta = (theta_lhs + theta_rhs + theta_delta) * 0.5;
            if (theta >= pi) {
                return generate_entire_sphere();
            }

            Real theta_rotate = theta - theta_lhs;
            Vector3f ra = cross(lhs.w, rhs.w);
            if (ra.length_squa() == 0) {
                return generate_entire_sphere();
            }
            return { Transform::generate_rotate(theta_rotate, ra)(lhs.w), std::cos(theta) };
        }

        bool is_empty() const noexcept {
            return std::isinf(cos_theta);
        }
        
        bool is_inside(const Vector3f &rhs) {
            return !is_empty() && dot(w, normalize(rhs)) >= cos_theta;
        }
    };
}
