#include <MatchOfflineRenderer/MatchOfflineRenderer.hpp>

using namespace MatchOfflineRenderer::math;

int main() {
    MatchOfflineRenderer::Initialize({
        // 指定是否计算 RGB转换表
        .generate_rgb_to_spectrum_table = false,
    });

    Vector3f v { 0, 1, 0 };
    Point3f p { 0, 2, 3 };
    Point3f p2 { 3, 2, 3 };

    MCH_INFO("v = {}\np = {}\np + v = {}\np - p2 = {}\n", v, p, p + v, p - p2)
    /*
    v = Vector3f { 0, 1, 0 }
    p = Point3f { 0, 2, 3 }
    p + v = Point3f { 0, 3, 3 }
    p - p2 = Vector3f { -3, 0, 0 }
    */

    Transform t = Transform::generate_translate({ 0, 2, 0 });
    MCH_INFO("p = {}, translate(p) = {}", p, t(p));
    /*
    p = Point3f { 0, 2, 3 }, translate(p) = Point3f { 0, 4, 3 }
    */

    // 弧度值
    Transform t_r = Transform::generate_rotate_x(radians(90));
    MCH_INFO("p = {}, rotate_x(p) = {}", p, t_r(p));
    /*
    p = Point3f { 0, 2, 3 }, rotate_x(p) = Point3f { 0, -3, 2 }
    */

    Transform t_v = Transform::generate_look_at({ 1, 1, 1 }, { 0, 0, 0 }, { 0, 1, 0 });
    MCH_INFO("p = {}, view(p) = {}", p, t_v(p));
    /*
    p = Point3f { 0, 2, 3 }, view(p) = Point3f { 2.12132, 0.408248, -1.1547 }
    */
    MCH_INFO("inv_view(p') = {}", t_v.apply_inverse(Point3f { 2.12132, 0.408248, -1.1547 }))
    /*
    inv_view(p') = Point3f { 1.19209e-07, 2, 3 }
    */

    Vector3f v2 = { 0, 0.4, 0.1};
    v2.normalize();
    MCH_INFO("normalize(v2) = {}", v2)
    auto theta = spherical_theta(v2);
    auto phi = spherical_phi(v2);
    MCH_INFO("theta = {}, phi = {}", theta, phi)
    // 弧度值
    /*
    theta = 1.3258177, phi = 1.5707964
    */
    MCH_INFO("v2 = {}", spherical_direction(sin(theta), cos(theta), phi))
    /*
    normalize(v2) = Vector3f { 0, 0.970142, 0.242536 }
    v2 = Vector3f { -4.24063e-08, 0.970142, 0.242536 }
    */

    // 点 - 点 = 向量
    // 点 + 向量 = 点
    // 向量 + 向量 = 向量
    // 向量 - 向量 = 向量
    // 点 + 点  不行

    // 采样波长, 会在360纳米到830纳米均匀的采样4个波长
    auto wavelengths =  MatchOfflineRenderer::spectra::SampledWavelengths::generate_sample_uniform(0.2);

    // 将标准D65光照转为sRGB色彩空间下的RGB值
    auto rgb = MatchOfflineRenderer::spectra::convert_spectrum_to_rgb(
        *MatchOfflineRenderer::spectra::piecewise_linear_cie_illum_d65,
        wavelengths,
        *MatchOfflineRenderer::spectra::colorspace_sRGB
    );
    MCH_INFO("R: {}, G: {}, B: {}", rgb.r, rgb.g, rgb.b);

    /*
    // 因为只有4个采样点,所以精度不高
    R: 1.3284488875245832, G: 0.8840898876193922, B: 2.270262453827748
    // 改为8000采样点, 很接近(1, 1, 1)了
    R: 0.9990772980420684, G: 0.9986944502601929, B: 1.0021058511109564
    */

    // 将sRGB色彩空间下的RGB(0.3, 0.5, 0.4)转换为反射光谱
    auto albedo = MatchOfflineRenderer::spectra::RGBAlbedoSpectrum(
        *MatchOfflineRenderer::spectra::colorspace_sRGB,
        { 0.3, 0.5, 0.4}
    );

    for (size_t lambda = 360; lambda < 830; lambda += 10) {
        MCH_INFO("在波长 {} 纳米下的反射率为 {}", lambda, albedo(lambda))
    }

    // 将反射光谱转换回sRGB色彩空间下的RGB值
    rgb = MatchOfflineRenderer::spectra::convert_spectrum_to_rgb(
        albedo, 
        wavelengths,
        *MatchOfflineRenderer::spectra::colorspace_sRGB
    );

    MCH_INFO("R: {}, G: {}, B: {}", rgb.r, rgb.g, rgb.b);
    /*
    (0.3, 0.5, 0.4)
    R: 0.37054229779340186, G: 0.4753282165780681, B: 0.3603767463851159
    改为4个采样点, 精度降低
    R: 0.5933764330246816, G: 0.440891496332491, B: 0.7718103597424845
    */

    // auto film = MatchOfflineRenderer::camera::RGBFilm {
    //     {},
    //     *MatchOfflineRenderer::spectra::colorspace_sRGB,
    //     1,
    //     false
    // };
    // film.resolution = { 1920, 1080 };
    // Medium medium;
    // auto c = MatchOfflineRenderer::camera::OrthographicCamera(
    //     {
    //         {},
    //         0,
    //         1,
    //         film,
    //         medium
    //     },
    //     { { 0, 0 }, { 1, 1 } },
    //     0, 1
    // );
    
    MatchOfflineRenderer::Destroy();
    return 0;
}
