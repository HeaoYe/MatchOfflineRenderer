#pragma once

#include <MatchOfflineRenderer/math/point.hpp>

namespace MatchOfflineRenderer::sample {
    struct Filter {
        math::Vector2f radius {};

        virtual math::Real integral() const noexcept = 0;

        virtual math::Real evaluate(const math::Point2f &point) const noexcept = 0;

        virtual ~Filter() noexcept = default;
    };
}
