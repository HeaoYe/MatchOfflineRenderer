#pragma once

#include <MatchOfflineRenderer/spectra/rgb_convert_types.hpp>
#include <MatchOfflineRenderer/spectra/color.hpp>
#include <MatchOfflineRenderer/core/loader.hpp>

namespace MatchOfflineRenderer::spectra {
    // 多项式光谱
    struct RGBSigmoidPolynomialSpectrum : public Spectrum {
        math::Real a {}, b {}, c {};

        RGBSigmoidPolynomialSpectrum() noexcept = default;

        RGBSigmoidPolynomialSpectrum(math::Real a, math::Real b, math::Real c) noexcept : a(a), b(b), c(c) {}

        math::Real operator()(math::Real lambda) const noexcept override { return math::sigmoid(math::evaluate_polynomial(lambda, c, b, a)); }

        math::Real max_value() const noexcept override {
            math::Real result = std::max((*this)(lambda_min), (*this)(lambda_max));
            math::Real lambda = -b / (2 * a);
            if (lambda > lambda_min && lambda < lambda_max) {
                result = std::max(result, (*this)(lambda));
            }
            return result;
        }
    };

    // RGB颜色与多项式系数的对应表,预计算得出
    struct RGBToSpectrumTable {
        const ZNodeArray<math::Real> z_nodes;
        const CoefficientArray<math::Real> coefficients;

        RGBToSpectrumTable(const ZNodeArray<math::Real> &z_nodes, const CoefficientArray<math::Real> &coefficients) : z_nodes(z_nodes), coefficients(coefficients) {}

        RGBSigmoidPolynomialSpectrum operator()(RGB rgb) const noexcept {
            if (rgb.r == rgb.g && rgb.g == rgb.b) {
                return { 0, 0, (rgb.r - 0.5) / std::sqrt(rgb.r * (1 - rgb.r)) };
            }

            math::Int maxc = rgb.max_component();
            math::Real z = rgb[maxc];
            math::Real x = rgb[(maxc + 1) % 3] * (rgb_to_spectrum_table_resolution - 1) / z;
            math::Real y = rgb[(maxc + 2) % 3] * (rgb_to_spectrum_table_resolution - 1) / z;

            math::Int xi = std::min<math::Int>(static_cast<math::Int>(x), rgb_to_spectrum_table_resolution - 2);
            math::Int yi = std::min<math::Int>(static_cast<math::Int>(y), rgb_to_spectrum_table_resolution - 2);
            math::Int zi = std::max<math::Int>(
                    std::upper_bound(z_nodes.begin(), z_nodes.end(), z) - z_nodes.begin() - 1,
                    0
                );
            math::Real dx = x - xi;
            math::Real dy = y - yi;
            math::Real dz = (z - z_nodes[zi]) / (z_nodes[zi + 1] - z_nodes[zi]);

            std::array<math::Real, 3> c;
            for (int i = 0; i < 3; ++i) {
                auto co = [&](int dx, int dy, int dz) {
                    return coefficients[maxc][zi + dz][yi + dy][xi + dx][i];
                };
                c[i] = math::lerp(
                    math::lerp(
                        math::lerp(co(0, 0, 0), co(1, 0, 0), dx),
                        math::lerp(co(0, 1, 0), co(1, 1, 0), dx),
                        dy
                    ),
                    math::lerp(
                        math::lerp(co(0, 0, 1), co(1, 0, 1), dx),
                        math::lerp(co(0, 1, 1), co(1, 1, 1), dx),
                        dy
                    ),
                    dz
                );
            }

            return { c[0], c[1], c[2] };
        }
    };
    
    extern std::unique_ptr<RGBToSpectrumTable> rgb_to_spectrum_table_sRGB;
    extern std::unique_ptr<RGBToSpectrumTable> rgb_to_spectrum_table_DCI_P3;
    extern std::unique_ptr<RGBToSpectrumTable> rgb_to_spectrum_table_Rec2020;
    extern std::unique_ptr<RGBToSpectrumTable> rgb_to_spectrum_table_ACES2065_1;

    // 预计算 RGB颜色与多项式系数的对应表 ,并写入文件中
    void GenerateRGBToSpectrumTables(const InitInfo &init_info);

    // 从预计算的RGB颜色与多项式系数的对应表加载数据
    void LoadRGBToSpectrumTables();
}
