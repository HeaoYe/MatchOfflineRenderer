#pragma once

#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/vector.hpp>

namespace MatchOfflineRenderer::math {
    // 光线
    
    struct Medium;

    struct RayPayload {
        Real time {};
        const Medium *medium { nullptr };
    };

    struct Ray {
        Point3f origin {};
        Vector3f direction {};
        RayPayload payload {};

        Ray() noexcept = default;
        
        Ray(const Point3f &origin, const Vector3f &direction, const RayPayload &payload = {}) noexcept : origin(origin), direction(direction), payload(payload) {}

        Point3f operator()(Real t) const noexcept { return origin + direction * t; }

        bool has_nan() const noexcept {
            return origin.has_nan() || direction.has_nan();
        }
    };

    struct RayDifferential : public Ray {
        bool has_differentials { false };
        Point3f origin_x {}, origin_y {};
        Vector3f direction_x {}, direction_y {};
        
        RayDifferential() noexcept = default;
        
        RayDifferential(const Point3f &origin, const Vector3f &direction, const RayPayload &payload = {}) noexcept : Ray(origin, direction, payload), has_differentials(false) {}
        
        RayDifferential(const Ray &ray) noexcept : Ray(ray), has_differentials(false) {}

        void scale_differentials(Real s) noexcept {
            origin_x = origin + (origin_x - origin) * s;
            origin_y = origin + (origin_y - origin) * s;
            direction_x = direction + (direction_x - direction) * s;
            direction_y = direction + (direction_y - direction) * s;
        }

        bool has_nan() const noexcept {
            return Ray::has_nan() || (
                has_differentials && 
                (origin_x.has_nan() || origin_y.has_nan()) || 
                (direction_x.has_nan() || direction_y.has_nan())
            );
        }
    };
}
