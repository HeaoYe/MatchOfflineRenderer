#pragma once

#include <MatchOfflineRenderer/spectra/spectrum.hpp>

namespace MatchOfflineRenderer::spectra {
    // 黑体辐射
    inline math::Real blackbody(math::Real lambda, math::Real t) {
        if (t <= 0) {
            return 0;
        }
        math::Real l = nm_to_m(1e-9);
        return (2 * h * c * c) / (std::pow<math::Real>(l, 5) * (std::exp((h * c) / (l * kb * t)) - 1));
    }
}
