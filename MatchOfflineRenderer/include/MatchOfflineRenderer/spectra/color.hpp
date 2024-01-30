#pragma once

#include <MatchOfflineRenderer/spectra/specified_spectrum.hpp>
#include <MatchOfflineRenderer/math/point.hpp>

namespace MatchOfflineRenderer::spectra {
    extern std::unique_ptr<spectra::DenselySampledSpectrum> densely_sampled_cie_x;
    extern std::unique_ptr<spectra::DenselySampledSpectrum> densely_sampled_cie_y;
    extern std::unique_ptr<spectra::DenselySampledSpectrum> densely_sampled_cie_z;

    // CIE公布的XYZ匹配函数
    inline const DenselySampledSpectrum &X() {
        return *densely_sampled_cie_x;
    }

    inline const DenselySampledSpectrum &Y() {
        return *densely_sampled_cie_y;
    }
    
    inline const DenselySampledSpectrum &Z() {
        return *densely_sampled_cie_z;
    }

    // XYZ颜色
    struct XYZ {
        math::Real x {}, y {}, z {};

        XYZ() noexcept = default;

        explicit XYZ(math::Real rhs) noexcept : x(rhs), y(rhs), z(rhs) {}

        XYZ(math::Real x, math::Real y, math::Real z) noexcept : x(x), y(y), z(z) {}

        explicit XYZ(const Spectrum &rhs) noexcept : x(inner_product(rhs, X()) / cie_y_integral), y(inner_product(rhs, Y()) / cie_y_integral), z(inner_product(rhs, Z()) / cie_y_integral) {}

        inline static XYZ from_xyy(math::Point2f xy, math::Real y = 1) {
            if (xy.y == 0) {
                return XYZ { 0 };
            }
            return { xy.x * y / xy.y, y, (1 - xy.x - xy.y) * y / xy.y };
        }

        math::Point2f xy() const {
            return { x / (x + y + z), y / (x + y + z) };
        }

        math::Real operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        math::Real &operator[](int i) noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        XYZ operator+(const XYZ &rhs) const noexcept {
            return { x + rhs.x, y + rhs.y, z + rhs.z };
        }

        XYZ operator-(const XYZ &rhs) const noexcept {
            return { x - rhs.x, y - rhs.y, z - rhs.z };
        }
        
        XYZ &operator+=(const XYZ &rhs) noexcept {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        XYZ &operator-=(const XYZ &rhs) noexcept {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        XYZ &operator-() noexcept {
            x = -x;
            y = -y;
            z = -z;
            return *this;
        }

        XYZ operator*(const XYZ &rhs) const noexcept {
            return { x * rhs.x, y * rhs.y, z * rhs.z };
        }

        XYZ operator/(const XYZ &rhs) const noexcept {
            return { x / rhs.x, y / rhs.y, z / rhs.z };
        }
        
        XYZ &operator*=(const XYZ &rhs) noexcept {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            return *this;
        }

        XYZ &operator/=(const XYZ &rhs) noexcept {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            return *this;
        }

        XYZ operator*(math::Real rhs) const noexcept {
            return { x * rhs, y * rhs, z * rhs };
        }

        XYZ operator/(math::Real rhs) const noexcept {
            return { x / rhs, y / rhs, z / rhs };
        }
        
        XYZ &operator*=(math::Real rhs) noexcept {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }

        XYZ &operator/=(math::Real rhs) noexcept {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }

        bool operator==(const XYZ &rhs) const noexcept { return x == rhs.x && y == rhs.y && z == rhs.z; }

        bool operator!=(const XYZ &rhs) const noexcept { return x != rhs.x || y != rhs.y || z != rhs.z; }
    };

    // 光谱转为XYZ颜色
    inline XYZ convert_sampled_spectrum_to_xyz(const SampledSpectrum &sampled_spectrum, const SampledWavelengths &wavelengths) noexcept {
        auto sampled_x = X().sample(wavelengths);
        auto sampled_y = Y().sample(wavelengths);
        auto sampled_z = Z().sample(wavelengths);

        auto sampled_pdf = wavelengths.pdf_as_sampled_spectrum();

        return XYZ {
            safe_div(sampled_x * sampled_spectrum, sampled_pdf).average(),
            safe_div(sampled_y * sampled_spectrum, sampled_pdf).average(),
            safe_div(sampled_z * sampled_spectrum, sampled_pdf).average()
        } / cie_y_integral;
    }

    inline XYZ convert_spectrum_to_xyz(const Spectrum &spectrum, const SampledWavelengths &wavelengths) noexcept {
        auto sampled_spectrum = spectrum.sample(wavelengths);
        return convert_sampled_spectrum_to_xyz(sampled_spectrum, wavelengths);
    }

    // RGB颜色
    struct RGB {
        math::Real r {}, g {}, b {};

        RGB() noexcept = default;

        explicit RGB(math::Real rhs) noexcept : r(rhs), g(rhs), b(rhs) {}

        RGB(math::Real r, math::Real g, math::Real b) noexcept : r(r), g(g), b(b) {}

        math::Real operator[](int i) const noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return r;
            if (i == 1) return g;
            return b;
        }

        math::Real &operator[](int i) noexcept {
            MCH_DASSERT((0 <= i) && (i < 3))
            if (i == 0) return r;
            if (i == 1) return g;
            return b;
        }

        math::Int max_component() const noexcept {
            return (r > g) ? ((r > b) ? 0 : 2) : ((g > b) ? 1 : 2);
        }

        math::Int min_component() const noexcept {
            return (r < g) ? ((r < b) ? 0 : 2) : ((g < b) ? 1 : 2);
        }

        RGB operator+(const RGB &rhs) const noexcept {
            return { r + rhs.r, g + rhs.g, b + rhs.b };
        }

        RGB operator-(const RGB &rhs) const noexcept {
            return { r - rhs.r, g - rhs.g, b - rhs.b };
        }
        
        RGB &operator+=(const RGB &rhs) noexcept {
            r += rhs.r;
            g += rhs.g;
            b += rhs.b;
            return *this;
        }

        RGB &operator-=(const RGB &rhs) noexcept {
            r -= rhs.r;
            g -= rhs.g;
            b -= rhs.b;
            return *this;
        }

        RGB &operator-() noexcept {
            r = -r;
            g = -g;
            b = -b;
            return *this;
        }

        RGB operator*(const RGB &rhs) const noexcept {
            return { r * rhs.r, g * rhs.g, b * rhs.b };
        }

        RGB operator/(const RGB &rhs) const noexcept {
            return { r / rhs.r, g / rhs.g, b / rhs.b };
        }
        
        RGB &operator*=(const RGB &rhs) noexcept {
            r *= rhs.r;
            g *= rhs.g;
            b *= rhs.b;
            return *this;
        }

        RGB &operator/=(const RGB &rhs) noexcept {
            r /= rhs.r;
            g /= rhs.g;
            b /= rhs.b;
            return *this;
        }

        RGB operator*(math::Real rhs) const noexcept {
            return { r * rhs, g * rhs, b * rhs };
        }

        RGB operator/(math::Real rhs) const noexcept {
            return { r / rhs, g / rhs, b / rhs };
        }
        
        RGB &operator*=(math::Real rhs) noexcept {
            r *= rhs;
            g *= rhs;
            b *= rhs;
            return *this;
        }

        RGB &operator/=(math::Real rhs) noexcept {
            r /= rhs;
            g /= rhs;
            b /= rhs;
            return *this;
        }

        bool operator==(const RGB &rhs) const noexcept { return r == rhs.r && g == rhs.g && b == rhs.b; }

        bool operator!=(const RGB &rhs) const noexcept { return r != rhs.r || g != rhs.g || b != rhs.b; }
    };

    inline RGB clamp_to_zero(const RGB &rhs) noexcept {
        return { std::max<math::Real>(rhs.r, 0), std::max<math::Real>(rhs.g, 0), std::max<math::Real>(rhs.b, 0) };
    }

    // 光谱转为RGB颜色涉及到了色彩空间,单独在别的头文件中实现
}
