#pragma once

#include <MatchOfflineRenderer/math/types.hpp>
#include <optional>

namespace MatchOfflineRenderer {
    // 正方形矩阵
    
    template <size_t N>
    struct SquareMatrix {
        Real m[N][N];

        SquareMatrix() noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = (i == j) ? 1 : 0;
                }
            }
        }
        
        SquareMatrix(const Real rhs[N][N]) noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = rhs[i][j];
                }
            }
        }
        
        SquareMatrix(const std::array<Real, N * N> &rhs) noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = rhs[i * N + j];
                }
            }
        }
        
        SquareMatrix(const std::array<std::array<Real, N>, N> &rhs) noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = rhs[i][j];
                }
            }
        }
        
        SquareMatrix(const std::vector<Real> &rhs) noexcept {
            MCH_DASSERT(rhs.size() == N * N)
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = rhs[i * N + j];
                }
            }
        }

        inline static SquareMatrix generate_zero_matrix() noexcept {
            SquareMatrix m;
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m.m[i][j] = 0;
                }
            }
            return m;
        }

        SquareMatrix operator+(const SquareMatrix &rhs) const noexcept {
            SquareMatrix result = *this;
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    result.m[i][j] += rhs.m[i][j];
                }
            }
            return result;
        }

        SquareMatrix operator*(Real rhs) const noexcept {
            SquareMatrix result = *this;
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    result.m[i][j] *= rhs;
                }
            }
            return result;
        }

        SquareMatrix operator/(Real rhs) const noexcept {
            rhs = Real { 1 } / rhs;
            MCH_DASSERT(!std::isnan(rhs))
            SquareMatrix result = *this;
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    result.m[i][j] *= rhs;
                }
            }
            return result;
        }

        bool operator==(const SquareMatrix &rhs) const noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    if (m[i][j] != rhs.m[i][j]) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool operator!=(const SquareMatrix &rhs) const noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    if (m[i][j] != rhs.m[i][j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        bool operator<(const SquareMatrix &rhs) const noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    if (m[i][j] < rhs.m[i][j]) return true;
                    if (m[i][j] > rhs.m[i][j]) return false;
                }
            }
            return false;
        }

        bool is_identify() const noexcept {
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    if (i == j) {
                        if (m[i][j] != 1) return false;
                    } else if (m[i][j] != 0) return false;
                }
            }
            return true;
        }

        const Real *operator[](int i) const noexcept {
            return m[i];
        }

        Real *operator[](int i) noexcept {
            return m[i];
        }

        template <typename R, typename T>
        R operator*(const T &rhs) const noexcept {
            R result;
            for (size_t i = 0; i < N; i ++) {
                result[i] = 0;
                for (size_t j = 0; j < N; j ++) {
                    result[i] += m[i][j] * rhs[j];
                }
            }
            return result;
        }
        
        Real determinant() const noexcept {
            Real minor12 = difference_of_products(m[1][1], m[2][2], m[1][2], m[2][1]);
            Real minor02 = difference_of_products(m[1][0], m[2][2], m[1][2], m[2][0]);
            Real minor01 = difference_of_products(m[1][0], m[2][1], m[1][1], m[2][0]);
            return std::fma<Real, Real, Real>(m[0][2], minor01, difference_of_products(m[0][0], minor12, m[0][1], minor02));
        }

        SquareMatrix &transpose() noexcept {
            SquareMatrix temp = *this;
            for (size_t i = 0; i < N; i ++) {
                for (size_t j = 0; j < N; j ++) {
                    m[i][j] = temp.m[j][i];
                }
            }
            return *this;
        }

        SquareMatrix &force_inverse() {
            auto inv = inverse(*this);
            MCH_ASSERT(inv.has_value())
            *this = inv.value();
            return *this;
        }
    };

    template <size_t N>
    inline SquareMatrix<N> operator*(Real lhs, const SquareMatrix<N> &rhs) noexcept {
        return rhs * lhs;
    }
   
    template <typename R, size_t N, typename T>
    inline R mul(const SquareMatrix<N> &lhs, const T &rhs) noexcept {
        return lhs * rhs;
    }

    template <size_t N, typename T>
    inline T operator*(const SquareMatrix<N> &lhs, const T &rhs) {
        return mul<T>(lhs, rhs);
    }

    template <>
    inline SquareMatrix<4> operator*(const SquareMatrix<4> &lhs, const SquareMatrix<4> &rhs) {
        SquareMatrix<4> result;
        for (size_t i = 0; i < 4; i ++) {
            for (size_t j = 0; j < 4; j ++) {
                result.m[i][j] = precise_inner_product(lhs.m[i][0], rhs.m[0][j], lhs.m[i][1], rhs.m[1][j], lhs.m[i][2], rhs.m[2][j], lhs.m[i][3], rhs.m[3][j]);
            }
        }
        return result;
    }

    template <>
    inline SquareMatrix<3> operator*(const SquareMatrix<3> &lhs, const SquareMatrix<3> &rhs) {
        SquareMatrix<3> result;
        for (size_t i = 0; i < 3; i ++) {
            for (size_t j = 0; j < 3; j ++) {
                result.m[i][j] = precise_inner_product(lhs.m[i][0], rhs.m[0][j], lhs.m[i][1], rhs.m[1][j], lhs.m[i][2], rhs.m[2][j]);
            }
        }
        return result;
    }

    template <size_t N>
    inline SquareMatrix<N> operator*(const SquareMatrix<N> &lhs, const SquareMatrix<N> &rhs) {
        SquareMatrix<N> result;
        for (size_t i = 0; i < N; i ++) {
            for (size_t j = 0; j < N; j ++) {
                result.m[i][j] = 0;
                for (size_t k = 0; k < N; k ++) {
                    result.m[i][j] = std::fma(lhs.m[i][k], rhs.m[k][j], result[i][j]);
                }
            }
        }
        return result;
    }
   
    template <size_t N>
    inline Real determinant(const SquareMatrix<N> &rhs) noexcept {
        return rhs.determinant();
    }
   
    template <size_t N>
    inline SquareMatrix<N> transpose(const SquareMatrix<N> &rhs) noexcept {
        SquareMatrix<N> result;
        for (size_t i = 0; i < N; i ++) {
            for (size_t j = 0; j < N; j ++) {
                result.m[i][j] = rhs.m[j][i];
            }
        }
        return result;
    }

    template <size_t N>
    inline std::optional<SquareMatrix<N>> inverse(const SquareMatrix<N> &rhs) noexcept {
        return {};
    }
   
    template <>
    inline std::optional<SquareMatrix<3>> inverse(const SquareMatrix<3> &rhs) noexcept {
        Real det = rhs.determinant();
        if (det == 0)
            return {};
        Real inv_det = Real { 1 } / det;

        SquareMatrix<3> result;

        result.m[0][0] = inv_det * difference_of_products(rhs.m[1][1], rhs.m[2][2], rhs.m[1][2], rhs.m[2][1]);
        result.m[1][0] = inv_det * difference_of_products(rhs.m[1][2], rhs.m[2][0], rhs.m[1][0], rhs.m[2][2]);
        result.m[2][0] = inv_det * difference_of_products(rhs.m[1][0], rhs.m[2][1], rhs.m[1][1], rhs.m[2][0]);
        result.m[0][1] = inv_det * difference_of_products(rhs.m[0][2], rhs.m[2][1], rhs.m[0][1], rhs.m[2][2]);
        result.m[1][1] = inv_det * difference_of_products(rhs.m[0][0], rhs.m[2][2], rhs.m[0][2], rhs.m[2][0]);
        result.m[2][1] = inv_det * difference_of_products(rhs.m[0][1], rhs.m[2][0], rhs.m[0][0], rhs.m[2][1]);
        result.m[0][2] = inv_det * difference_of_products(rhs.m[0][1], rhs.m[1][2], rhs.m[0][2], rhs.m[1][1]);
        result.m[1][2] = inv_det * difference_of_products(rhs.m[0][2], rhs.m[1][0], rhs.m[0][0], rhs.m[1][2]);
        result.m[2][2] = inv_det * difference_of_products(rhs.m[0][0], rhs.m[1][1], rhs.m[0][1], rhs.m[1][0]);

        return result;
    }

    template <>
    inline std::optional<SquareMatrix<4>> inverse(const SquareMatrix<4> &rhs) noexcept {
        // Via: https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
        // (c) Google, Apache license.

        // For 4x4 do not compute the adjugate as the transpose of the cofactor
        // matrix, because this results in extra work. Several calculations can be
        // shared across the sub-determinants.
        //
        // This approach is explained in David Eberly's Geometric Tools book,
        // excerpted here:
        //   http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
        Real s0 = difference_of_products(rhs.m[0][0], rhs.m[1][1], rhs.m[1][0], rhs.m[0][1]);
        Real s1 = difference_of_products(rhs.m[0][0], rhs.m[1][2], rhs.m[1][0], rhs.m[0][2]);
        Real s2 = difference_of_products(rhs.m[0][0], rhs.m[1][3], rhs.m[1][0], rhs.m[0][3]);

        Real s3 = difference_of_products(rhs.m[0][1], rhs.m[1][2], rhs.m[1][1], rhs.m[0][2]);
        Real s4 = difference_of_products(rhs.m[0][1], rhs.m[1][3], rhs.m[1][1], rhs.m[0][3]);
        Real s5 = difference_of_products(rhs.m[0][2], rhs.m[1][3], rhs.m[1][2], rhs.m[0][3]);

        Real c0 = difference_of_products(rhs.m[2][0], rhs.m[3][1], rhs.m[3][0], rhs.m[2][1]);
        Real c1 = difference_of_products(rhs.m[2][0], rhs.m[3][2], rhs.m[3][0], rhs.m[2][2]);
        Real c2 = difference_of_products(rhs.m[2][0], rhs.m[3][3], rhs.m[3][0], rhs.m[2][3]);

        Real c3 = difference_of_products(rhs.m[2][1], rhs.m[3][2], rhs.m[3][1], rhs.m[2][2]);
        Real c4 = difference_of_products(rhs.m[2][1], rhs.m[3][3], rhs.m[3][1], rhs.m[2][3]);
        Real c5 = difference_of_products(rhs.m[2][2], rhs.m[3][3], rhs.m[3][2], rhs.m[2][3]);

        Real det = precise_inner_product(s0, c5, -s1, c4, s2, c3, s3, c2, s5, c0, -s4, c1);
        if (det == 0)
            return {};
        Real inv_det = Real { 1 } / det;

        Real inv[4][4] = {
            {
                inv_det * precise_inner_product(rhs.m[1][1], c5, rhs.m[1][3], c3, -rhs.m[1][2], c4),
                inv_det * precise_inner_product(-rhs.m[0][1], c5, rhs.m[0][2], c4, -rhs.m[0][3], c3),
                inv_det * precise_inner_product(rhs.m[3][1], s5, rhs.m[3][3], s3, -rhs.m[3][2], s4),
                inv_det * precise_inner_product(-rhs.m[2][1], s5, rhs.m[2][2], s4, -rhs.m[2][3], s3)
            }, {
                inv_det * precise_inner_product(-rhs.m[1][0], c5, rhs.m[1][2], c2, -rhs.m[1][3], c1),
                inv_det * precise_inner_product(rhs.m[0][0], c5, rhs.m[0][3], c1, -rhs.m[0][2], c2),
                inv_det * precise_inner_product(-rhs.m[3][0], s5, rhs.m[3][2], s2, -rhs.m[3][3], s1),
                inv_det * precise_inner_product(rhs.m[2][0], s5, rhs.m[2][3], s1, -rhs.m[2][2], s2)
            }, {
                inv_det * precise_inner_product(rhs.m[1][0], c4, rhs.m[1][3], c0, -rhs.m[1][1], c2),
                inv_det * precise_inner_product(-rhs.m[0][0], c4, rhs.m[0][1], c2, -rhs.m[0][3], c0),
                inv_det * precise_inner_product(rhs.m[3][0], s4, rhs.m[3][3], s0, -rhs.m[3][1], s2),
                inv_det * precise_inner_product(-rhs.m[2][0], s4, rhs.m[2][1], s2, -rhs.m[2][3], s0)
            }, {
                inv_det * precise_inner_product(-rhs.m[1][0], c3, rhs.m[1][1], c1, -rhs.m[1][2], c0),
                inv_det * precise_inner_product(rhs.m[0][0], c3, rhs.m[0][2], c0, -rhs.m[0][1], c1),
                inv_det * precise_inner_product(-rhs.m[3][0], s3, rhs.m[3][1], s1, -rhs.m[3][2], s0),
                inv_det * precise_inner_product(rhs.m[2][0], s3, rhs.m[2][2], s0, -rhs.m[2][1], s1)
            }
        };

        return SquareMatrix<4>(inv);
    }
   
    template <size_t N>
    inline SquareMatrix<N> force_inverse(const SquareMatrix<N> &rhs) {
        auto inv = inverse(rhs);
        MCH_ASSERT(inv.has_value())
        return inv.value();
    }

    template <size_t N>
    std::optional<SquareMatrix<N>> linear_least_squares(const Real A[][N], const Real B[][N], size_t rows) {
        SquareMatrix<N> AtA = SquareMatrix<N>::generate_zero_matrix();
        SquareMatrix<N> AtB = SquareMatrix<N>::generate_zero_matrix();

        for (size_t i = 0; i < N; i ++) {
            for (size_t j = 0; j < N; j ++) {
                for (size_t r = 0; r < rows; r ++) {
                    AtA.m[i][j] += A[r][i] * A[r][j];
                    AtB.m[i][j] += A[r][i] * B[r][j];
                }
            }
        }

        auto AtAi = inverse(AtA);
        if (!AtAi.has_value())
            return {};
        return transpose(AtAi.value() * AtB);
    }
}
