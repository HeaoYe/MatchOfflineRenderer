#pragma once

#include <MatchOfflineRenderer/math/types.hpp>

namespace MatchOfflineRenderer::spectra {
    constexpr int rgb_to_spectrum_table_resolution = 64;

    template <typename T>
    using ZNodeArray = std::array<T, rgb_to_spectrum_table_resolution>;
    template <typename T>
    using CoefficientArray = std::array<std::array<std::array<std::array<std::array<T, 3>, rgb_to_spectrum_table_resolution>, rgb_to_spectrum_table_resolution>, rgb_to_spectrum_table_resolution>, 3>;

    extern const ZNodeArray<math::Real> rgb_to_spectrum_table_sRGB_z_nodes;
    extern const CoefficientArray<math::Real> rgb_to_spectrum_table_sRGB_coefficients;
    extern const ZNodeArray<math::Real> rgb_to_spectrum_table_DCI_P3_z_nodes;
    extern const CoefficientArray<math::Real> rgb_to_spectrum_table_DCI_P3_coefficients;
    extern const ZNodeArray<math::Real> rgb_to_spectrum_table_Rec2020_z_nodes;
    extern const CoefficientArray<math::Real> rgb_to_spectrum_table_Rec2020_coefficients;
    extern const ZNodeArray<math::Real> rgb_to_spectrum_table_ACES2065_1_z_nodes;
    extern const CoefficientArray<math::Real> rgb_to_spectrum_table_ACES2065_1_coefficients;
}
