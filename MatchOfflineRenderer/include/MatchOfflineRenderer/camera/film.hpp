#pragma once

#include <MatchOfflineRenderer/camera/image.hpp>
#include <MatchOfflineRenderer/math/vector.hpp>
#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/normal.hpp>
#include <MatchOfflineRenderer/math/bounds.hpp>
#include <MatchOfflineRenderer/math/interaction.hpp>
#include <MatchOfflineRenderer/spectra/spectrum.hpp>
#include <MatchOfflineRenderer/spectra/colorspace.hpp>
#include <MatchOfflineRenderer/sample/filter.hpp>
#include <atomic>

namespace MatchOfflineRenderer::camera {
    template <typename T>
    inline T project_reflectance(const spectra::Spectrum &reflectance, const spectra::Spectrum &illum, const spectra::Spectrum &x, const spectra::Spectrum &y, const spectra::Spectrum &z) noexcept {
        T result { 0, 0, 0 };
        math::Real y_integral { 0 };
        for (math::Real lambda = spectra::lambda_min; lambda <= spectra::lambda_max; lambda ++) {
            y_integral += y(lambda) * illum(lambda);
            result[0] += x(lambda) * illum(lambda) * reflectance(lambda);
            result[1] += y(lambda) * illum(lambda) * reflectance(lambda);
            result[2] += z(lambda) * illum(lambda) * reflectance(lambda);
        }
        return result;
    }

    inline math::SquareMatrix<3> white_balance(math::Point2f src_white, math::Point2f dst_white) noexcept {
        auto src_lms = spectra::xyz_to_lms(spectra::XYZ::from_xyy(src_white));
        auto dst_lms = spectra::xyz_to_lms(spectra::XYZ::from_xyy(dst_white));
        
        auto lms_correct = math::SquareMatrix<3>::generate_diag_matrix({
            dst_lms.x / src_lms.x,
            dst_lms.y / src_lms.y,
            dst_lms.z / src_lms.z
        });
        return spectra::lms_to_xyz_matrix * lms_correct * spectra::xyz_to_lms_matrix;
    }

    // 像素传感器
    struct PixelSensor {
        spectra::DenselySampledSpectrum r;
        spectra::DenselySampledSpectrum g;
        spectra::DenselySampledSpectrum b;
        math::Real imaging_ratio {};
        math::SquareMatrix<3> sensor_rgb_to_xyz {};

        static constexpr math::Int swatch_reflectance_count = 24;
        static spectra::PiecewiseLinearSpectrum swatch_reflectances[swatch_reflectance_count];

        PixelSensor(const spectra::Spectrum &r, const spectra::Spectrum &g, const spectra::Spectrum &b, const spectra::RGBColorSpace &output_colorspace, math::Real imaging_ratio, const spectra::Spectrum *sensor_illum = nullptr) : r(r), g(g), b(b), imaging_ratio(imaging_ratio) {
            math::Real rgb_camera[swatch_reflectance_count][3];
            math::Real xyz_output[swatch_reflectance_count][3];
            math::Real sensor_white_g = spectra::inner_product(*sensor_illum, g);
            math::Real sensor_white_y = spectra::inner_product(*sensor_illum, spectra::Y());

            for (size_t i = 0; i < swatch_reflectance_count; i ++) {
                auto rgb = project_reflectance<spectra::RGB>(swatch_reflectances[i], *sensor_illum, r, g, b);
                auto xyz = project_reflectance<spectra::XYZ>(swatch_reflectances[i], *sensor_illum, spectra::X(), spectra::Y(), spectra::Z());
                for (size_t c = 0; c < 3; c ++) {
                    rgb_camera[i][c] = rgb[c];
                    xyz_output[i][c] = xyz[c] * sensor_white_y / sensor_white_g;
                }
            }

            auto m = math::linear_least_squares(rgb_camera, xyz_output, swatch_reflectance_count);
            if (!m.has_value()) {
                MCH_FATAL("Error {} {}", __FILE__, __LINE__)
            }
            sensor_rgb_to_xyz = m.value();
        }

        PixelSensor(const spectra::RGBColorSpace &output_colorspace, math::Real imaging_ratio, const spectra::Spectrum *sensor_illum = nullptr) noexcept : r(spectra::X()), g(spectra::Y()), b(spectra::Z()), imaging_ratio(imaging_ratio) {
            if (sensor_illum != nullptr) {
                spectra::XYZ src_white { *sensor_illum };
                spectra::XYZ dst_white { output_colorspace.illuminant };
                sensor_rgb_to_xyz = white_balance(src_white.xy(), dst_white.xy());
            }
        }

        spectra::RGB convert_to_sensor_rgb(const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &sampled_wavelengths) const noexcept {
            auto sampled_spectrum_div_pdf = spectra::safe_div(radiance, sampled_wavelengths.pdf_as_sampled_spectrum());
            return {
                (r.sample(sampled_wavelengths) * sampled_spectrum_div_pdf).average(),
                (g.sample(sampled_wavelengths) * sampled_spectrum_div_pdf).average(),
                (b.sample(sampled_wavelengths) * sampled_spectrum_div_pdf).average()
            };
        }
    };

    // 光线命中的表面信息
    struct VisibleSurface {
        math::Point3f point {};
        math::Normal3f normal {}, shading_normal {};
        math::Point2f uv {};
        math::Real time { 0 };
        math::Vector3f dpdx {}, dpdy {};
        spectra::SampledSpectrum albedo { 0 };

        bool initialized { false };

        VisibleSurface() noexcept = default;

        VisibleSurface(const math::SurfaceInteraction &surface_interaction, const spectra::SampledSpectrum &albedo, const spectra::SampledWavelengths &wavelengths) noexcept;

        operator bool() const noexcept { return initialized; }
    };

    struct FilmCreateInfo {
        math::Point2i resolution {};
        math::Bounds2i pixel_bounds {};
        math::Real diagonal {};
        const PixelSensor *sensor { nullptr };
        const sample::Filter *filter { nullptr };
    };

    // 相机胶片
    struct Film {
        math::Point2i resolution {};
        math::Bounds2i pixel_bounds {};
        math::Real diagonal {};
        std::unique_ptr<const PixelSensor> sensor {};
        std::unique_ptr<const sample::Filter> filter {};

        Film(const FilmCreateInfo &info) noexcept : resolution(info.resolution), pixel_bounds(info.pixel_bounds), diagonal(info.diagonal * 0.001), sensor(info.sensor), filter(info.filter) {}

        virtual math::Bounds2f sample_bounds() const noexcept {
            math::Vector2f radius = filter->radius;
            return {
                pixel_bounds.min - radius + math::Vector2f(0.5f, 0.5f),
                pixel_bounds.max + radius - math::Vector2f(0.5f, 0.5f)
            };
        }

        virtual spectra::SampledWavelengths sample_wavelengths(math::Real t) const noexcept {
            return spectra::SampledWavelengths::generate_sample_visiable(t);
        }

        virtual bool uses_visible_surface() const noexcept = 0;
        virtual void add_sample(const math::Point2i &film_point, const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths, const VisibleSurface *visibleSurface, math::Real weight) noexcept = 0;
        virtual void add_splat(const math::Point2f &point, const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths) noexcept = 0;
        virtual spectra::RGB get_pixel_rgb(const math::Point2i &point, math::Real splat_scale = 1) const noexcept = 0;
        virtual spectra::RGB convert_to_output_rgb(const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths) const noexcept = 0;
        // virtual void write_image(const ImageMetadata &metadata, math::Real splat_scale = 1) noexcept = 0;
    };

    // RGB胶片
    struct RGBFilm : public Film {
        struct RGBPixel {
            double rgb_sum[3] = {0., 0., 0.};
            double weight_sum = 0.;
            std::atomic<uint64_t> rgb_splat[3];
        };

        const spectra::RGBColorSpace &colorspace;
        math::Real max_component_value { std::numeric_limits<math::Real>::infinity() };
        math::Real filter_integral {};
        bool write_fp16 {};
        math::SquareMatrix<3> sensor_rgb_to_output_rgb;
        math::Int pixel_width {};
        std::vector<RGBPixel> pixels {};

        RGBFilm(const FilmCreateInfo &info, const spectra::RGBColorSpace &colorspace, math::Real max_component_value, bool write_fp16) noexcept : Film(info), colorspace(colorspace), pixels(pixel_bounds.area()), max_component_value(max_component_value), write_fp16(write_fp16) {
            pixel_width = pixel_bounds.diagonal().x;
            filter_integral = filter->integral();
            sensor_rgb_to_output_rgb = colorspace.xyz_to_rgb * sensor->sensor_rgb_to_xyz;
        }

        virtual bool uses_visible_surface() const noexcept override { return false; }

        virtual void add_sample(const math::Point2i &film_point, const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths, const VisibleSurface *visibleSurface, math::Real weight) noexcept override {
            auto rgb = sensor->convert_to_sensor_rgb(radiance, wavelengths);
            auto max_value = rgb[rgb.max_component()];
            if (max_value > max_component_value) {
                rgb *= max_component_value / max_value;
            }
            auto &pixel = pixels[film_point.y * pixel_width + film_point.x];
            for (size_t c = 0; c < 3; c ++) {
                pixel.rgb_sum[c] += weight * rgb[c];
            }
            pixel.weight_sum += weight;
        }

        virtual void add_splat(const math::Point2f &point, const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths) noexcept override {
            auto rgb = sensor->convert_to_sensor_rgb(radiance, wavelengths);
            auto max_value = rgb[rgb.max_component()];
            if (max_value > max_component_value) {
                rgb *= max_component_value / max_value;
            }
            auto discrete_point = point + math::Vector2f { 0.5, 0.5 };
            auto radius = filter->radius;
            math::Bounds2i splat_bounds {
                math::Point2i { math::floor(discrete_point - radius) },
                math::Point2i { math::floor(discrete_point + radius) } + math::Vector2i { 1, 1 }
            };
            splat_bounds = math::Bounds2i::generate_intersect(splat_bounds, pixel_bounds);
            for (auto splat_point : splat_bounds) {
                math::Real wt = filter->evaluate(math::Point2f { point - splat_point - math::Vector2f { 0.5, 0.5 } });
                if (wt != 0) {
                    auto &pixel = pixels[splat_point.y * pixel_width + splat_point.x];
                    for (size_t c = 0; c < 3; c ++) {
                        uint64_t old_bits = pixel.rgb_splat[c];
                        uint64_t new_bits;
                        do {
                            new_bits = math::double_to_bits(math::bits_to_double(old_bits) + double(wt * rgb[c]));
                        } while (!std::atomic_compare_exchange_weak(&pixel.rgb_splat[c], &old_bits, new_bits));
                    }
                }
            }
        }

        virtual spectra::RGB get_pixel_rgb(const math::Point2i &point, math::Real splat_scale = 1) const noexcept override {
            auto &pixel = pixels[point.y * pixel_width + point.x];
            double r = pixel.rgb_sum[0];
            double g = pixel.rgb_sum[1];
            double b = pixel.rgb_sum[2];
            if (pixel.weight_sum != 0) {
                r /= pixel.weight_sum;
                g /= pixel.weight_sum;
                b /= pixel.weight_sum;
            }
            r += double(splat_scale) * math::bits_to_double(pixel.rgb_splat[0]) / filter_integral;
            g += double(splat_scale) * math::bits_to_double(pixel.rgb_splat[1]) / filter_integral;
            b += double(splat_scale) * math::bits_to_double(pixel.rgb_splat[2]) / filter_integral;

            return sensor_rgb_to_output_rgb * spectra::RGB { r, g, b };
        }

        virtual spectra::RGB convert_to_output_rgb(const spectra::SampledSpectrum &radiance, const spectra::SampledWavelengths &wavelengths) const noexcept override {
            auto sensor_rgb = sensor->convert_to_sensor_rgb(radiance, wavelengths);
            return sensor_rgb_to_output_rgb * sensor_rgb;
        }
    };
}
