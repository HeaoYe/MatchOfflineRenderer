#pragma once

#include <MatchOfflineRenderer/math/types.hpp>

namespace MatchOfflineRenderer::spectra {
    // 光谱相关常量
    constexpr math::Real c = math::Real { 299792458.0 };
    constexpr math::Real h = math::Real { 6.62606957e-34 };
    constexpr math::Real kb = math::Real { 1.3806488e-23 };

    constexpr math::Real lambda_min = 360;
    constexpr math::Real lambda_max = 830;

    // 采样数
    // 4个采样数是兼顾效率与质量的最优选择
    constexpr math::Int spectrum_sample_count = 4;

    constexpr math::Real cie_y_integral = 106.856895;

    inline math::Real nm_to_m(math::Real rhs) {
        return rhs * math::Real { 1e-9 };
    }
}
