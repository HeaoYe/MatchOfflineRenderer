#include <MatchOfflineRenderer/spectra/loader.hpp>
#include <MatchOfflineRenderer/spectra/cie_constants.hpp>
#include <MatchOfflineRenderer/spectra/color.hpp>
#include <MatchOfflineRenderer/spectra/colorspace.hpp>

namespace MatchOfflineRenderer::spectra {
    std::unique_ptr<DenselySampledSpectrum> densely_sampled_cie_x { nullptr };
    std::unique_ptr<DenselySampledSpectrum> densely_sampled_cie_y { nullptr };
    std::unique_ptr<DenselySampledSpectrum> densely_sampled_cie_z { nullptr };

    std::unique_ptr<PiecewiseLinearSpectrum> piecewise_linear_cie_illum_d65 { nullptr };
    std::unique_ptr<PiecewiseLinearSpectrum> piecewise_linear_cie_aces_illum_d60 { nullptr };

    std::unique_ptr<RGBColorSpace> colorspace_sRGB { nullptr };
    std::unique_ptr<RGBColorSpace> colorspace_DCI_P3 { nullptr };
    std::unique_ptr<RGBColorSpace> colorspace_Rec2020 { nullptr };
    std::unique_ptr<RGBColorSpace> colorspace_ACES2065_1 { nullptr };

    void Initialize(const InitInfo &init_info) {
        // 加载CIE公布的测量数据

        // XYZ匹配函数
        PiecewiseLinearSpectrum linear_x { spectra::cie_lambdas, spectra::cie_x };
        PiecewiseLinearSpectrum linear_y { spectra::cie_lambdas, spectra::cie_y };
        PiecewiseLinearSpectrum linear_z { spectra::cie_lambdas, spectra::cie_z };

        densely_sampled_cie_x = std::make_unique<DenselySampledSpectrum>(linear_x);
        densely_sampled_cie_y = std::make_unique<DenselySampledSpectrum>(linear_y);
        densely_sampled_cie_z = std::make_unique<DenselySampledSpectrum>(linear_z);

        // 标准光照光谱
        piecewise_linear_cie_illum_d65 = std::make_unique<PiecewiseLinearSpectrum>(convert_from_cie_data(cie_illum_d65, true));
        piecewise_linear_cie_aces_illum_d60 = std::make_unique<PiecewiseLinearSpectrum>(convert_from_cie_data(cie_aces_illum_d60, true));
        
        // 色彩空间
        colorspace_sRGB = std::make_unique<RGBColorSpace>(
            math::Point2f(.64, .33), math::Point2f(.3, .6), math::Point2f(.15, .06),
            *piecewise_linear_cie_illum_d65
        );
        
        colorspace_DCI_P3 = std::make_unique<RGBColorSpace>(
            math::Point2f(.68, .32), math::Point2f(.265, .690), math::Point2f(.15, .06),
            *piecewise_linear_cie_illum_d65
        );
        
        colorspace_Rec2020 = std::make_unique<RGBColorSpace>(
            math::Point2f(.708, .292), math::Point2f(.170, .797), math::Point2f(.131, .046),
            *piecewise_linear_cie_illum_d65
        );

        colorspace_ACES2065_1 = std::make_unique<RGBColorSpace>(
            math::Point2f(.7347, .2653), math::Point2f(0., 1.), math::Point2f(.0001, -.077),
            *piecewise_linear_cie_aces_illum_d60
        );

        // RGB转换表
        if (init_info.generate_rgb_to_spectrum_table) {
            GenerateRGBToSpectrumTables(init_info);
        } else {
            LoadRGBToSpectrumTables();
        }

        // 绑定色彩空间与转换表
        colorspace_sRGB->setup_rgb_to_spectrum_table(rgb_to_spectrum_table_sRGB.get());
        colorspace_DCI_P3->setup_rgb_to_spectrum_table(rgb_to_spectrum_table_DCI_P3.get());
        colorspace_Rec2020->setup_rgb_to_spectrum_table(rgb_to_spectrum_table_Rec2020.get());
        colorspace_ACES2065_1->setup_rgb_to_spectrum_table(rgb_to_spectrum_table_ACES2065_1.get());
    }

    void Destroy() {
        colorspace_ACES2065_1.reset();
        colorspace_Rec2020.reset();
        colorspace_DCI_P3.reset();
        colorspace_sRGB.reset();

        piecewise_linear_cie_aces_illum_d60.reset();
        piecewise_linear_cie_illum_d65.reset();

        densely_sampled_cie_z.reset();
        densely_sampled_cie_y.reset();
        densely_sampled_cie_x.reset();
    }
}
