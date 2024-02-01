#pragma once

#include <MatchOfflineRenderer/spectra/spectrum.hpp>
#include <MatchOfflineRenderer/spectra/blackbody.hpp>

namespace MatchOfflineRenderer::spectra {
    // 常量光谱
    struct ConstantSpectrum : public Spectrum {
        math::Real c {};

        ConstantSpectrum(math::Real c) noexcept : c(c) {}

        math::Real operator()(math::Real lambda) const noexcept override { return c; }

        math::Real max_value() const noexcept override { return -std::numeric_limits<math::Real>::infinity(); }
    };

    // Densely光谱,从lambda_min到lambda_max,每隔1nm采样一个波长对应的光强度
    struct DenselySampledSpectrum : public Spectrum {
        math::Int lambda_min {};
        math::Int lambda_max {};
        std::vector<math::Real> values {};

        DenselySampledSpectrum() noexcept = default;

        DenselySampledSpectrum(const Spectrum &spectrum, int lambda_min = ::MatchOfflineRenderer::spectra::lambda_min, int lambda_max = ::MatchOfflineRenderer::spectra::lambda_max) noexcept : lambda_min(lambda_min), lambda_max(lambda_max), values(lambda_max - lambda_min + 1) {
            for (int lambda = lambda_min; lambda <= lambda_max; lambda ++) {
                values[lambda - lambda_min] = spectrum(lambda);
            }
        }

        math::Real operator()(math::Real lambda) const noexcept override {
            math::Int offset = std::lround<math::Int>(lambda) - lambda_min;
            if (offset < 0 || offset >= values.size()) {
                return 0;
            }
            return values[offset];
        }

        math::Real max_value() const noexcept override { return *std::max_element(values.begin(), values.end()); }
    };

    // PiecewiseLinear光谱,由多个(lambda, value)组成的连续光谱
    struct PiecewiseLinearSpectrum : public Spectrum {
        std::vector<math::Real> lambdas {};
        std::vector<math::Real> values {};

        PiecewiseLinearSpectrum(const Spectrum &spectrum, const std::vector<math::Real> &lambdas) noexcept : PiecewiseLinearSpectrum(lambdas, ([&]() {
            std::vector<math::Real> values(lambdas.size());
            for (size_t i = 0; i < lambdas.size(); i ++) {
                values[i] = spectrum(lambdas[i]);
            }
            return values;
        })()) {}

        PiecewiseLinearSpectrum(const std::vector<math::Real> &lambdas, const std::vector<math::Real> &values) noexcept {
            MCH_ASSERT(lambdas.size() == values.size())

            std::vector<std::pair<math::Real, math::Real>> lambdas_values(lambdas.size());
            for (size_t i = 0; i < lambdas.size(); i ++) {
                lambdas_values[i] = std::make_pair(lambdas[i], values[i]);
            }
            std::sort(lambdas_values.begin(), lambdas_values.end(), [](auto &lhs, auto &rhs) {
                return lhs.first < rhs.first;
            });

            this->lambdas.reserve(lambdas_values.size());
            this->values.reserve(lambdas_values.size());
            for (size_t i = 0; i < lambdas_values.size(); i ++) {
                this->lambdas.push_back(lambdas_values[i].first);
                this->values.push_back(lambdas_values[i].second);
            }
        }

        template <size_t N>
        PiecewiseLinearSpectrum(const std::array<math::Real, N> &lambdas, const std::array<math::Real, N> &values) noexcept {
            std::array<std::pair<math::Real, math::Real>, N> lambdas_values;
            for (size_t i = 0; i < lambdas.size(); i ++) {
                lambdas_values[i] = std::make_pair(lambdas[i], values[i]);
            }
            std::sort(lambdas_values.begin(), lambdas_values.end(), [](auto &lhs, auto &rhs) {
                return lhs.first < rhs.first;
            });

            this->lambdas.reserve(lambdas_values.size());
            this->values.reserve(lambdas_values.size());
            for (size_t i = 0; i < lambdas_values.size(); i ++) {
                this->lambdas.push_back(lambdas_values[i].first);
                this->values.push_back(lambdas_values[i].second);
            }
        }

        math::Real operator()(math::Real lambda) const noexcept override {
            if (lambdas.empty() || lambda < lambdas.front() || lambda > lambdas.back()) {
                return 0;
            }

            for (size_t i = 1; i < lambdas.size(); i ++) {
                if ((lambdas[i - 1] <= lambda) && (lambda < lambdas[i])) {
                    math::Real t = (lambda - lambdas[i - 1]) / (lambdas[i] - lambdas[i - 1]);
                    return (1 - t) * values[i - 1] + t * values[i];
                }
            }

            if (lambda == lambdas.back()) {
                return values.back();
            }

            MCH_ASSERT(false)
            return 0;
        }

        math::Real max_value() const noexcept override { return *std::max_element(values.begin(), values.end()); }
    };

    // 黑体辐射光谱
    struct BlackbodySpectrum : public Spectrum {
        math::Real t {};
        math::Real normalization_factor {};

        BlackbodySpectrum(math::Real t) noexcept : t(t) {
            math::Real lambda_max = math::Real { 2.8977721e-3 } / t;
            normalization_factor = math::Real { 1 } / blackbody(nm_to_m(lambda_max), t);
        };

        math::Real operator()(math::Real lambda) const noexcept override {
            return blackbody(lambda, t) * normalization_factor;
        }

        math::Real max_value() const noexcept override { return 1; }
    };

    // CIE公布的D65光谱和D60光谱
    extern std::unique_ptr<spectra::PiecewiseLinearSpectrum> piecewise_linear_cie_illum_d65;
    extern std::unique_ptr<spectra::PiecewiseLinearSpectrum> piecewise_linear_cie_aces_illum_d60;

    template <size_t N>
    PiecewiseLinearSpectrum convert_from_interleaved_data(const std::array<math::Real, N> &samples, bool normalize) noexcept {
        extern std::unique_ptr<spectra::DenselySampledSpectrum> densely_sampled_cie_y;
        MCH_DASSERT(N % 2 == 0)

        constexpr size_t n = N / 2;
        std::vector<math::Real> lambdas;
        std::vector<math::Real> values;
        lambdas.reserve(n + 2);
        values.reserve(n + 2);

        if (samples[0] > lambda_min) {
            lambdas.push_back(lambda_min - 1);
            values.push_back(samples[1]);
        }
        for (size_t i = 0; i < n; i ++) {
            lambdas.push_back(samples[2 * i]);
            values.push_back(samples[2 * i + 1]);
            MCH_DASSERT(i == 0 || lambdas.back() > lambdas[lambdas.size() - 2]);
        }
        if (lambdas.back() < lambda_max) {
            lambdas.push_back(lambda_max + 1);
            values.push_back(values.back());
        }

        PiecewiseLinearSpectrum spectrum { lambdas, values };

        if (normalize) {
            math::Real scale = cie_y_integral / inner_product(spectrum, *densely_sampled_cie_y);
            for (auto &value : spectrum.values) {
                value *= scale;
            }
        }

        return spectrum;
    }
}
