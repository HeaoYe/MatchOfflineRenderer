#pragma once

#include <MatchOfflineRenderer/math/point.hpp>

namespace MatchOfflineRenderer::sample {
    inline math::Point2f sample_uniform_disk_polar(const math::Point2f &t) {
        math::Real r = std::sqrt(t[0]);
        math::Real theta = 2 * math::pi * t[1];
        return { r * std::cos(theta), r * std::sin(theta) };
    }

    inline math::Point2f sample_uniform_disk_concentric(const math::Point2f &t) {
        auto offset = t * 2 - math::Vector2f { 1, 1 };
        if (offset.x == 0 && offset.y == 0) {
            return { 0, 0 };
        }
        math::Real theta, r;
        if (std::abs(offset.x) > std::abs(offset.y)) {
            r = offset.x;
            theta = math::pi_div4 * (offset.y / offset.x);
        } else {
            r = offset.y;
            theta = math::pi_div2 - math::pi_div4 * (offset.x / offset.y);
        }
        return math::Point2f { std::cos(theta), std::sin(theta) } * r;
    }
}
