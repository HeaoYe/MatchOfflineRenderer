#pragma once

#include <MatchOfflineRenderer/spectra/colorspace.hpp>

namespace MatchOfflineRenderer::spectra {
    // 反射RGB光谱,每个波长的值在[0, 1]的范围内, 遵循能量守恒, 且波长对应的值连续
    struct RGBAlbedoSpectrum : public Spectrum {
        RGBSigmoidPolynomialSpectrum sigmoid_polynomial_spectrum {};

        RGBAlbedoSpectrum(const RGBColorSpace &colorspace, const RGB &rgb) {
            MCH_DASSERT(
                rgb.r >= 0 && rgb.r <= 1 &&
                rgb.g >= 0 && rgb.g <= 1 &&
                rgb.b >= 0 && rgb.b <= 1
            )
            sigmoid_polynomial_spectrum = colorspace.convert_to_sigmoid_polynomial_spectrum(rgb);
        }

        math::Real operator()(math::Real lambda) const noexcept override { return sigmoid_polynomial_spectrum(lambda); }

        math::Real max_value() const noexcept override { return sigmoid_polynomial_spectrum.max_value(); }
    };

    // 无边界RGB光谱,每个波长的值可以超过1, 且波长对应的值连续
    struct RGBUnboundedSpectrum : public RGBAlbedoSpectrum {
        math::Real scale { 1 };

        RGBUnboundedSpectrum(const RGBColorSpace &colorspace, const RGB &rgb) noexcept : 
            RGBAlbedoSpectrum(
                colorspace,
                [&]() {
                    math::Real m = std::max({ rgb.r, rgb.g, rgb.b });
                    scale = 2 * m;
                    return (scale == 0) ? RGB { 0, 0, 0 } : (rgb / scale);
                }()
            ) {}

        math::Real operator()(math::Real lambda) const noexcept override { return scale * sigmoid_polynomial_spectrum(lambda); }

        math::Real max_value() const noexcept override { return scale * sigmoid_polynomial_spectrum.max_value(); }
    };

    // 光照RGB光谱,每个波长的值可以超过1, 且波长对应的值连续
    struct RGBIlluminantSpectrum : public RGBUnboundedSpectrum {
        const DenselySampledSpectrum &illuminant;

        RGBIlluminantSpectrum(const RGBColorSpace &colorspace, const RGB &rgb) noexcept : RGBUnboundedSpectrum(colorspace, rgb), illuminant(colorspace.illuminant) {}

        math::Real operator()(math::Real lambda) const noexcept override {
            return scale * sigmoid_polynomial_spectrum(lambda) * illuminant(lambda);
        }
    };
}
