#pragma once

#include <MatchOfflineRenderer/spectra/constants.hpp>

namespace MatchOfflineRenderer::spectra {
    // 对光谱的采样
    struct SampledSpectrum {
        std::array<math::Real, spectrum_sample_count> values {};

        explicit SampledSpectrum(math::Real c) noexcept { values.fill(c); }

        SampledSpectrum(const std::array<math::Real, spectrum_sample_count> &values) noexcept : values(values) {}

        bool is_zero() const noexcept {
            for (auto value : values) {
                if (value != 0) {
                    return false;
                }
            }
            return true;
        }

        SampledSpectrum operator+(const SampledSpectrum &rhs) const noexcept {
            SampledSpectrum result(values);
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                result.values[i] += rhs.values[i];
            }
            return result;
        }

        SampledSpectrum operator-(const SampledSpectrum &rhs) const noexcept {
            SampledSpectrum result(values);
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                result.values[i] -= rhs.values[i];
            }
            return result;
        }

        SampledSpectrum &operator+=(const SampledSpectrum &rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++)
                values[i] += rhs.values[i];
            return *this;
        }

        SampledSpectrum &operator-=(const SampledSpectrum &rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] -= rhs.values[i];
            }
            return *this;
        }

        SampledSpectrum &operator-() noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] = -values[i];
            }
            return *this;
        }

        SampledSpectrum operator*(const SampledSpectrum &rhs) const noexcept {
            SampledSpectrum result(values);
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                result.values[i] *= rhs.values[i];
            }
            return result;
        }

        SampledSpectrum operator/(const SampledSpectrum &rhs) const noexcept {
            SampledSpectrum result(values);
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                result.values[i] /= rhs.values[i];
            }
            return result;
        }

        SampledSpectrum &operator*=(const SampledSpectrum &rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++)
                values[i] *= rhs.values[i];
            return *this;
        }

        SampledSpectrum &operator/=(const SampledSpectrum &rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] /= rhs.values[i];
            }
            return *this;
        }

        SampledSpectrum &safe_div(const SampledSpectrum &rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] = (rhs.values[i] == 0) ? (values[i] / rhs.values[i]) : 0;
            }
            return *this;
        }

        SampledSpectrum &sqrt() noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] = std::sqrt(values[i]);
            }
            return *this;
        }

        SampledSpectrum &clamp(math::Real min, math::Real max) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] = std::clamp<math::Real>(values[i], min, max);
            }
            return *this;
        }

        SampledSpectrum &pow(math::Real rhs) noexcept {
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                values[i] = std::pow<math::Real>(values[i], rhs);
            }
            return *this;
        }

        math::Real min_component() const noexcept {
            return *std::min_element(values.begin(), values.end());
        }

        math::Real max_component() const noexcept {
            return *std::max_element(values.begin(), values.end());
        }

        math::Real average() const noexcept {
            math::Real sum = 0;
            for (auto value : values) {
                sum += value;
            }
            return sum / spectrum_sample_count;
        }
    };

    inline SampledSpectrum safe_div(const SampledSpectrum &lhs, const SampledSpectrum &rhs) noexcept {
        SampledSpectrum result(0);
        for (size_t i = 0; i < spectrum_sample_count; i ++) {
            result.values[i] = (rhs.values[i] == 0) ? 0 : (lhs.values[i] / rhs.values[i]);
        }
        return result;
    }

    inline SampledSpectrum lerp(const SampledSpectrum &lhs, const SampledSpectrum &rhs, math::Real t) noexcept {
        SampledSpectrum result(0);
        for (size_t i = 0; i < spectrum_sample_count; i ++) {
            result.values[i] = math::lerp(lhs.values[i], rhs.values[i], t);
        }
        return result;
    }

    inline SampledSpectrum sqrt(const SampledSpectrum &rhs) noexcept {
        SampledSpectrum result(0);
        for (size_t i = 0; i < spectrum_sample_count; i ++) {
            result.values[i] = std::sqrt(rhs.values[i]);
        }
        return result;
    }

    inline SampledSpectrum clamp(const SampledSpectrum &rhs, math::Real min, math::Real max) noexcept {
        SampledSpectrum result(0);
        for (size_t i = 0; i < spectrum_sample_count; i ++) {
            result.values[i] = std::clamp<math::Real>(rhs.values[i], min, max);
        }
        return result;
    }

    inline SampledSpectrum pow(const SampledSpectrum &lhs, math::Real rhs) noexcept {
        SampledSpectrum result(0);
        for (size_t i = 0; i < spectrum_sample_count; i ++) {
            result.values[i] = std::pow<math::Real>(lhs.values[i], rhs);
        }
        return result;
    }

    // 对波长的采样
    struct SampledWavelengths {
        std::array<math::Real, spectrum_sample_count> lambda {};
        std::array<math::Real, spectrum_sample_count> pdf {};

        SampledSpectrum pdf_as_sampled_spectrum() const noexcept { return { pdf }; }

        bool is_terminated() const noexcept {
            for (size_t i = 1; i < spectrum_sample_count; i ++) {
                if (pdf[i] != 0) {
                    return false;
                }
            }
            return true;
        }

        void terminate() noexcept {
            if (is_terminated()) {
                return;
            }
            for (size_t i = 1; i < spectrum_sample_count; i ++) {
                pdf[i] = 0;
            }
            pdf[0] /= spectrum_sample_count;
        }

        inline static SampledWavelengths generate_sample_uniform(math::Real t, math::Real lambda_min = ::MatchOfflineRenderer::spectra::lambda_min, math::Real lambda_max = ::MatchOfflineRenderer::spectra::lambda_max) noexcept {
            SampledWavelengths result {};

            result.lambda[0] = math::lerp(lambda_min, lambda_max, t);

            math::Real delta = (lambda_max - lambda_min) / spectrum_sample_count;
            for (size_t i = 1; i < spectrum_sample_count; i ++) {
                result.lambda[i] = result.lambda[i - 1] + delta;
                if (result.lambda[i] > lambda_max) {
                    result.lambda[i] = lambda_min + (result.lambda[i] - lambda_max);
                }
            }

            result.pdf.fill(1 / (lambda_max - lambda_min));

            return result;
        }
    };

    // 光谱基类
    struct Spectrum {
        virtual math::Real operator()(math::Real lambda) const noexcept = 0;
        virtual math::Real max_value() const noexcept = 0;
        virtual ~Spectrum() = default;

        SampledSpectrum sample(const SampledWavelengths &wavelengths) const noexcept {
            SampledSpectrum result { 0 };
            for (size_t i = 0; i < spectrum_sample_count; i ++) {
                result.values[i] = (*this)(wavelengths.lambda[i]);
            }
            return result;
        }
    };

    inline math::Real inner_product(const Spectrum &lhs, const Spectrum &rhs) noexcept {
        math::Real integral = 0;
        for (math::Real lambda = lambda_min; lambda <= lambda_max; lambda ++) {
            auto a = lhs(lambda), b = rhs(lambda);
            integral += lhs(lambda) * rhs(lambda);
        }
        return integral;
    }
}
