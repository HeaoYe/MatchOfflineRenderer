#include <MatchOfflineRenderer/MatchOfflineRenderer.hpp>
#include <MatchOfflineRenderer/math/tuple.hpp>

using namespace MatchOfflineRenderer;

int main() {
    Initialize();

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
    
    Destroy();
    return 0;
}
