#pragma once

#include <MatchOfflineRenderer/math/bounds.hpp>
#include <MatchOfflineRenderer/math/sphere.hpp>

namespace MatchOfflineRenderer::shape {
    struct Shape {
        math::Bounds3f bounds {};
        math::DirectionCone normal_bounds {};
    };
}
