#pragma once

#include <MatchOfflineRenderer/math/matrix.hpp>
#include <MatchOfflineRenderer/math/vector.hpp>
#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/normal.hpp>
#include <MatchOfflineRenderer/math/ray.hpp>
#include <MatchOfflineRenderer/math/bounds.hpp>

namespace MatchOfflineRenderer::math {
    // 坐标系

    struct CoordinateSystem {
        Vector3f normal { 0, 1, 0 };
        Vector3f tangent { 1, 0, 0 };
        Vector3f bitangent { 0, 0, 1 };

        CoordinateSystem() noexcept = default;

        CoordinateSystem(const Vector3f &normal, const Vector3f &tangent, const Vector3f &bitangent) noexcept : normal(normal), tangent(tangent), bitangent(bitangent) {}

        Vector3f to_local(const Vector3f &rhs) const noexcept {
            return { dot(rhs, tangent), dot(rhs, normal), dot(rhs, bitangent) };
        }

        Normal3f ToLocal(const Normal3f &rhs) const noexcept {
            return Normal3f { dot(rhs, tangent), dot(rhs, normal), dot(rhs, bitangent) };
        }

        Vector3f from_local(const Vector3f &rhs) const noexcept {
            return tangent * rhs.x + normal * rhs.y + bitangent * rhs.z;
        }
        
        Normal3f from_local(const Normal3f &rhs) const noexcept {
            return Normal3f { tangent * rhs.x + normal * rhs.y + bitangent * rhs.z };
        }

        inline static CoordinateSystem from_xy(const Vector3f &x, const Vector3f &y) noexcept {
            return CoordinateSystem(y, x, cross(x, y));
        }

        inline static CoordinateSystem from_xz(const Vector3f &x, const Vector3f &z) noexcept {
            return CoordinateSystem(cross(z, x), x, z);
        }

        inline static CoordinateSystem from_yz(const Vector3f &y, const Vector3f &z) noexcept {
            return CoordinateSystem(y, cross(y, z), y);
        }

        inline static CoordinateSystem from_x(const Vector3f &x) noexcept {
            Real sign = std::copysign(Real { 1 }, x.z);
            Real a = -1 / (sign + x.z);
            Real b = x.x * x.y * a;
            return {
                { 1 + sign * x.x * x.x * a, sign * b, -sign * x.x },
                x,
                { b, sign + x.y * x.y * a, -x.y },
            };
        }

        inline static CoordinateSystem from_x(const Normal3f &x) noexcept {
            Real sign = std::copysign(Real { 1 }, x.z);
            Real a = -1 / (sign + x.z);
            Real b = x.x * x.y * a;
            return {
                { 1 + sign * x.x * x.x * a, sign * b, -sign * x.x },
                Vector3f(x),
                { b, sign + x.y * x.y * a, -x.y },
            };
        }

        inline static CoordinateSystem from_y(const Vector3f &y) noexcept {
            Real sign = std::copysign(Real { 1 }, y.z);
            Real a = -1 / (sign + y.z);
            Real b = y.x * y.y * a;
            return {
                y,
                { 1 + sign * y.x * y.x * a, sign * b, -sign * y.x },
                { b, sign + y.y * y.y * a, -y.y },
            };
        }

        inline static CoordinateSystem from_y(const Normal3f &y) noexcept {
            Real sign = std::copysign(Real { 1 }, y.z);
            Real a = -1 / (sign + y.z);
            Real b = y.x * y.y * a;
            return {
                Vector3f(y),
                { 1 + sign * y.x * y.x * a, sign * b, -sign * y.x },
                { b, sign + y.y * y.y * a, -y.y },
            };
        }

        inline static CoordinateSystem from_z(const Vector3f &z) noexcept {
            Real sign = std::copysign(Real { 1 }, z.z);
            Real a = -1 / (sign + z.z);
            Real b = z.x * z.y * a;
            return {
                { b, sign + z.y * z.y * a, -z.y },
                { 1 + sign * z.x * z.x * a, sign * b, -sign * z.x },
                z,
            };
        }

        inline static CoordinateSystem from_z(const Normal3f &z) noexcept {
            Real sign = std::copysign(Real { 1 }, z.z);
            Real a = -1 / (sign + z.z);
            Real b = z.x * z.y * a;
            return {
                { b, sign + z.y * z.y * a, -z.y },
                { 1 + sign * z.x * z.x * a, sign * b, -sign * z.x },
                Vector3f(z),
            };
        }
    };

    // 变换
    struct Transform {
        SquareMatrix<4> matrix {}, inv_matrix {};

        Transform() noexcept = default;

        Transform(const SquareMatrix<4> &rhs) noexcept : matrix(rhs) {
            auto inv = ::MatchOfflineRenderer::math::inverse(rhs);
            if (inv.has_value()) {
                inv_matrix = inv.value();
            } else {
                Real nan = std::numeric_limits<Real>::has_signaling_NaN ? std::numeric_limits<Real>::signaling_NaN() : std::numeric_limits<Real>::quiet_NaN();
                for (size_t i = 0; i < 4; i ++) {
                    for (size_t j = 0; j < 4; j ++) {
                        inv_matrix.m[i][j] = nan;
                    }
                }
            }
        }

        Transform(const Real rhs[4][4]) noexcept : Transform(SquareMatrix<4> { rhs }) {}

        Transform(const SquareMatrix<4> &matrix, const SquareMatrix<4> &inv_matrix) noexcept : matrix(matrix), inv_matrix(inv_matrix) {}

        explicit Transform(const CoordinateSystem &coord) noexcept : 
            Transform(std::array<Real, 16> {
                coord.tangent.x, coord.tangent.y, coord.tangent.z, 0,
                coord.normal.x, coord.normal.y, coord.normal.z, 0,
                coord.bitangent.x, coord.bitangent.y, coord.bitangent.z, 0,
                0, 0, 0, 1
            }) {}
 
        Transform &inverse() noexcept {
            std::swap(matrix, inv_matrix);
            return *this;
        }

        Transform &transpose() noexcept {
            matrix.transpose();
            inv_matrix.transpose();
            return *this;
        }

        bool operator==(const Transform &rhs) const noexcept { return rhs.matrix == matrix; }
        
        bool operator!=(const Transform &rhs) const noexcept { return rhs.matrix != matrix; }
        
        bool is_identify() const noexcept { return matrix.is_identify(); }
        
        bool has_scale(Real tolerance = 1e-3f) const noexcept {
            Real la2 = ((*this)(Vector3f { 1, 0, 0 })).length_squa();
            Real lb2 = ((*this)(Vector3f { 0, 1, 0 })).length_squa();
            Real lc2 = ((*this)(Vector3f { 0, 0, 1 })).length_squa();
            return (std::abs(la2 - 1) > tolerance ||
                    std::abs(lb2 - 1) > tolerance ||
                    std::abs(lc2 - 1) > tolerance);
        }

        Transform operator*(const Transform &rhs) const noexcept {
            return { matrix * rhs.matrix, rhs.inv_matrix * inv_matrix };
        }

        bool swaps_handedness() const noexcept {
            SquareMatrix<3> s {
                std::array<Real, 9> {
                    matrix.m[0][0], matrix.m[0][1], matrix.m[0][2],
                    matrix.m[1][0], matrix.m[1][1], matrix.m[1][2],
                    matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]
                }
            };
            return s.determinant() < 0;
        }

        template <typename T>
        Point3<T> operator()(const Point3<T> &rhs) const noexcept {
            T xp = matrix.m[0][0] * rhs.x + matrix.m[0][1] * rhs.y + matrix.m[0][2] * rhs.z + matrix.m[0][3];
            T yp = matrix.m[1][0] * rhs.x + matrix.m[1][1] * rhs.y + matrix.m[1][2] * rhs.z + matrix.m[1][3];
            T zp = matrix.m[2][0] * rhs.x + matrix.m[2][1] * rhs.y + matrix.m[2][2] * rhs.z + matrix.m[2][3];
            T wp = matrix.m[3][0] * rhs.x + matrix.m[3][1] * rhs.y + matrix.m[3][2] * rhs.z + matrix.m[3][3];
            if (wp == 1) return { xp, yp, zp };
            return Point3<T> { xp, yp, zp } / wp;
        }

        template <typename T>
        Point3<T> apply_inverse(const Point3<T> &rhs) const noexcept {
            T xp = inv_matrix.m[0][0] * rhs.x + inv_matrix.m[0][1] * rhs.y + inv_matrix.m[0][2] * rhs.z + inv_matrix.m[0][3];
            T yp = inv_matrix.m[1][0] * rhs.x + inv_matrix.m[1][1] * rhs.y + inv_matrix.m[1][2] * rhs.z + inv_matrix.m[1][3];
            T zp = inv_matrix.m[2][0] * rhs.x + inv_matrix.m[2][1] * rhs.y + inv_matrix.m[2][2] * rhs.z + inv_matrix.m[2][3];
            T wp = inv_matrix.m[3][0] * rhs.x + inv_matrix.m[3][1] * rhs.y + inv_matrix.m[3][2] * rhs.z + inv_matrix.m[3][3];
            if (wp == 1) return { xp, yp, zp };
            return Point3<T> { xp, yp, zp } / wp;
        }

        template <typename T>
        Vector3<T> operator()(const Vector3<T> &rhs) const noexcept {
            return {
                matrix.m[0][0] * rhs.x + matrix.m[0][1] * rhs.y + matrix.m[0][2] * rhs.z,
                matrix.m[1][0] * rhs.x + matrix.m[1][1] * rhs.y + matrix.m[1][2] * rhs.z,
                matrix.m[2][0] * rhs.x + matrix.m[2][1] * rhs.y + matrix.m[2][2] * rhs.z
            };
        }

        template <typename T>
        Vector3<T> apply_inverse(const Vector3<T> &rhs) const noexcept {
            return {
                inv_matrix.m[0][0] * rhs.x + inv_matrix.m[0][1] * rhs.y + inv_matrix.m[0][2] * rhs.z,
                inv_matrix.m[1][0] * rhs.x + inv_matrix.m[1][1] * rhs.y + inv_matrix.m[1][2] * rhs.z,
                inv_matrix.m[2][0] * rhs.x + inv_matrix.m[2][1] * rhs.y + inv_matrix.m[2][2] * rhs.z
            };
        }

        template <typename T>
        Normal3<T> operator()(const Normal3<T> &rhs) const noexcept {
            return {
                inv_matrix.m[0][0] * rhs.x + inv_matrix.m[1][0] * rhs.y + inv_matrix.m[2][0] * rhs.z,
                inv_matrix.m[0][1] * rhs.x + inv_matrix.m[1][1] * rhs.y + inv_matrix.m[2][1] * rhs.z,
                inv_matrix.m[0][2] * rhs.x + inv_matrix.m[1][2] * rhs.y + inv_matrix.m[2][2] * rhs.z
            };
        }

        Ray operator()(const Ray &rhs, Real *tMax) const noexcept {
            Point3fi origin = (*this)(Point3fi(rhs.origin));
            Vector3f direction = (*this)(rhs.direction);
            if (Real length_squa = direction.length_squa(); length_squa > 0) {
                // Real dt = Dot(abs(direction), origin.Error()) / length_squa;
                // origin += direction * dt;
                // if (tMax)
                //     *tMax -= dt;
            }
            return { Point3f(origin), direction, rhs.payload };
        }
 
        Bounds3f operator()(const Bounds3f &rhs) const noexcept {
            Bounds3f result;
            for (int i = 0; i < 8; ++i)
                result = Bounds3f::generate_union(result, (*this)(rhs.corner(i)));
            return result;
        }

        inline static Transform generate_translate(const Vector3f &delta) noexcept {
            return {
                std::array<Real, 16> {
                    1, 0, 0, delta.x,
                    0, 1, 0, delta.y,
                    0, 0, 1, delta.z,
                    0, 0, 0, 1
                },
                std::array<Real, 16> {
                    1, 0, 0, -delta.x,
                    0, 1, 0, -delta.y,
                    0, 0, 1, -delta.z,
                    0, 0, 0, 1
                }
            };
        }

        inline static Transform generate_scale(Real x, Real y, Real z) noexcept {
            return {
                std::array<Real, 16> {
                    x, 0, 0, 0,
                    0, y, 0, 0,
                    0, 0, z, 0,
                    0, 0, 0, 1
                },
                std::array<Real, 16> {
                    1 / x, 0, 0, 0,
                    0, 1 / y, 0, 0,
                    0, 0, 1 / z, 0,
                    0, 0, 0, 1
                }
            };
        }

        inline static Transform generate_rotate_x(Real theta) noexcept {
            Real sin_theta = std::sin(theta);
            Real cos_theta = std::cos(theta);
            return {
                std::array<Real, 16> {
                    1, 0, 0, 0,
                    0, cos_theta, -sin_theta, 0,
                    0, sin_theta, cos_theta, 0,
                    0, 0, 0, 1
                }
            };
        }

        inline static Transform generate_rotate_y(Real theta) noexcept {
            Real sin_theta = std::sin(theta);
            Real cos_theta = std::cos(theta);
            return {
                std::array<Real, 16> {
                    cos_theta, 0, sin_theta, 0,
                    0, 1, 0, 0,
                    -sin_theta, 0, cos_theta, 0,
                    0, 0, 0, 1
                }
            };
        }

        inline static Transform generate_rotate_z(Real theta) noexcept {
            Real sin_theta = std::sin(theta);
            Real cos_theta = std::cos(theta);
            return {
                std::array<Real, 16> {
                    cos_theta, -sin_theta, 0, 0,
                    sin_theta, cos_theta, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1
                }
            };
        }

        inline static Transform generate_rotate(Real sin_theta, Real cos_theta, const Vector3f &axis) noexcept {
            auto a = normalize(axis);
            SquareMatrix<4> m;
            
            m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cos_theta;
            m.m[0][1] = a.x * a.y * (1 - cos_theta) - a.z * sin_theta;
            m.m[0][2] = a.x * a.z * (1 - cos_theta) + a.y * sin_theta;
            m.m[0][3] = 0;

            m.m[1][0] = a.x * a.y * (1 - cos_theta) + a.z * sin_theta;
            m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cos_theta;
            m.m[1][2] = a.y * a.z * (1 - cos_theta) - a.x * sin_theta;
            m.m[1][3] = 0;
            
            m.m[2][0] = a.x * a.z * (1 - cos_theta) - a.y * sin_theta;
            m.m[2][1] = a.y * a.z * (1 - cos_theta) + a.x * sin_theta;
            m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cos_theta;
            m.m[2][3] = 0;

            return { m };
        }
        
        inline static Transform generate_rotate(Real theta, const Vector3f &axis) noexcept {
            return generate_rotate(std::sin(theta), std::cos(theta), axis);
        }
        
        inline static Transform generate_rotate_from_to(const Vector3f &from, const Vector3f &to) noexcept {
            Vector3f refl;
            if (std::abs(from.x) < 0.72f && std::abs(to.x) < 0.72f) {
                refl = Vector3f(1, 0, 0);
            } else if (std::abs(from.y) < 0.72f && std::abs(to.y) < 0.72f) {
                refl = Vector3f(0, 1, 0);
            } else {
                refl = Vector3f(0, 0, 1);
            }

            Vector3f u = refl - from, v = refl - to;
            SquareMatrix<4> result;
            for (size_t i = 0; i < 3; i ++) {
                for (size_t j = 0; j < 3; j ++) {
                    result.m[i][j] = ((i == j) ? 1 : 0) - 
                        2 / dot(u, u) * u[i] * u[j] -
                        2 / dot(v, v) * v[i] * v[j] +
                        4 * dot(u, v) / (dot(u, u) * dot(v, v)) * v[i] * u[j];
                }
            }
            return { result };
        }

        inline static Transform generate_look_at(const Point3f &pos, const Point3f &look, const Vector3f &up) noexcept {
            SquareMatrix<4> camera_to_world;

            camera_to_world.m[0][3] = pos.x;
            camera_to_world.m[1][3] = pos.y;
            camera_to_world.m[2][3] = pos.z;
            camera_to_world.m[3][3] = 1;

            Vector3f dir = normalize(look - pos);
            Vector3f right = normalize(cross(normalize(up), dir));
            Vector3f new_up = cross(dir, right);

            camera_to_world.m[0][0] = right.x;
            camera_to_world.m[1][0] = right.y;
            camera_to_world.m[2][0] = right.z;
            camera_to_world.m[3][0] = 0.;
            camera_to_world.m[0][1] = new_up.x;
            camera_to_world.m[1][1] = new_up.y;
            camera_to_world.m[2][1] = new_up.z;
            camera_to_world.m[3][1] = 0.;
            camera_to_world.m[0][2] = dir.x;
            camera_to_world.m[1][2] = dir.y;
            camera_to_world.m[2][2] = dir.z;
            camera_to_world.m[3][2] = 0.;

            SquareMatrix<4> world_to_camera = force_inverse(camera_to_world);
            return { world_to_camera, camera_to_world };
        }
    };

    inline Transform inverse(const Transform &rhs) noexcept {
        return { rhs.inv_matrix, rhs.matrix };
    }
    
    inline Transform transpose(const Transform &rhs) {
        return { transpose(rhs.matrix), transpose(rhs.inv_matrix) };
    }
}
