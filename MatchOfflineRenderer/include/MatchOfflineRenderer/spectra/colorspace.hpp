#pragma once

#include <MatchOfflineRenderer/spectra/rgb_convert.hpp>
#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/matrix.hpp>

namespace MatchOfflineRenderer::spectra {
    // 色彩空间
    struct RGBColorSpace {
        math::Point2f r {}, g {}, b {}, w {};
        DenselySampledSpectrum illuminant {};
        const RGBToSpectrumTable *rgb_to_spectrum_table { nullptr };
        math::SquareMatrix<3> xyz_to_rgb {}, rgb_to_xyz {};

        RGBColorSpace(math::Point2f r, math::Point2f g, math::Point2f b, Spectrum &illuminant) noexcept : r(r), g(g), b(b), illuminant(illuminant) {
            XYZ w_xyz { illuminant };
            w = w_xyz.xy();

            XYZ r_xyz = XYZ::from_xyy(r);
            XYZ g_xyz = XYZ::from_xyy(g);
            XYZ b_xyz = XYZ::from_xyy(b);

            math::SquareMatrix<3> rgb = std::array {
                r_xyz.x, g_xyz.x, b_xyz.x,
                r_xyz.y, g_xyz.y, b_xyz.y,
                r_xyz.z, g_xyz.z, b_xyz.z
            };

            XYZ c = math::force_inverse(rgb) * w_xyz;
            rgb_to_xyz = rgb * math::SquareMatrix<3>::generate_diag_matrix({ c.x, c.y, c.z });
            xyz_to_rgb = math::force_inverse(rgb_to_xyz);
        }

        void setup_rgb_to_spectrum_table(const RGBToSpectrumTable *rgb_to_spectrum_table) {
            this->rgb_to_spectrum_table = rgb_to_spectrum_table;
        }
        
        // XYZ颜色与RGB颜色的互相转换

        RGB convert_to_rgb(const XYZ &xyz) const noexcept { return math::mul<RGB>(xyz_to_rgb, xyz); }
        
        XYZ convert_to_xyz(const RGB &rgb) const noexcept { return math::mul<XYZ>(rgb_to_xyz, rgb); }

        // RGB颜色转化为光谱
        RGBSigmoidPolynomialSpectrum convert_to_sigmoid_polynomial_spectrum(const RGB &rgb) const noexcept {
            MCH_DASSERT(rgb_to_spectrum_table != nullptr)
            return (*rgb_to_spectrum_table)(clamp_to_zero(rgb));
        }

        bool operator==(const RGBColorSpace &rhs) const noexcept {
            return (r == rhs.r && g == rhs.g && b == rhs.b && w == rhs.w &&
                    rgb_to_spectrum_table == rhs.rgb_to_spectrum_table);
        }

        bool operator!=(const RGBColorSpace &rhs) const noexcept {
            return (r != rhs.r || g != rhs.g || b != rhs.b || w != rhs.w ||
                    rgb_to_spectrum_table != rhs.rgb_to_spectrum_table);
        }
    };

    inline math::SquareMatrix<3> convert_rgb_colorspace(const RGBColorSpace &from, const RGBColorSpace &to) noexcept {
        if (from == to) {
            return {};
        }
        return to.xyz_to_rgb * from.rgb_to_xyz;
    }

    // 光谱转RGB颜色

    inline RGB convert_sampled_spectrum_to_rgb(const SampledSpectrum &sampled_spectrum, const SampledWavelengths &wavelengths, const RGBColorSpace &rgb_colorspace) noexcept {
        XYZ xyz = convert_sampled_spectrum_to_xyz(sampled_spectrum, wavelengths);
        return rgb_colorspace.convert_to_rgb(xyz);
    }

    inline RGB convert_spectrum_to_rgb(const Spectrum &spectrum, const SampledWavelengths &wavelengths, const RGBColorSpace &rgb_colorspace) noexcept {
        return convert_sampled_spectrum_to_rgb(spectrum.sample(wavelengths), wavelengths, rgb_colorspace);
    }

    // 目前支持的色彩空间
    enum class RGBColorSpaceType {
        esRGB,
        eDCI_P3,
        eRec2020,
        eACES2065_1,
    };

    extern std::unique_ptr<RGBColorSpace> colorspace_sRGB;
    extern std::unique_ptr<RGBColorSpace> colorspace_DCI_P3;
    extern std::unique_ptr<RGBColorSpace> colorspace_Rec2020;
    extern std::unique_ptr<RGBColorSpace> colorspace_ACES2065_1;
}
