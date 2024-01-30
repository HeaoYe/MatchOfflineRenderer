#pragma once

#include <MatchOfflineRenderer/math/interval.hpp>
#include <MatchOfflineRenderer/math/vector.hpp>
#include <MatchOfflineRenderer/math/point.hpp>
#include <MatchOfflineRenderer/math/normal.hpp>

namespace MatchOfflineRenderer::math {
    struct SurfaceInteraction;
    struct MediumInteraction;
    struct MediumInterface;

    struct Medium {};
    struct PhaseFunction {};

    // 交互信息
    struct Interaction {
        Point3fi pointi {};
        Vector3f wo {};
        Real time {};

        Normal3f normal {};
        Point2f uv {};

        const MediumInterface *mediumInterface { nullptr };
        Medium medium {};

        Interaction(const Point3fi &pointi, const Normal3f &normal, const Point2f &uv, const Vector3f &wo, Real time)
            : pointi(pointi), normal(normal), uv(uv), wo(wo), time(time) {}

        Interaction(Point3f point, Vector3f wo, Real time, Medium medium)
            : pointi(point), wo(wo), time(time), medium(medium) {}

        Point3f point() const noexcept { return Point3f(pointi); }
        
        bool is_surface_interaction() const noexcept { return normal != Normal3f(0, 0, 0); }

        bool is_medium_interaction() const noexcept { return !is_surface_interaction(); }
        
        const SurfaceInteraction &as_surface() const noexcept {
            MCH_DASSERT(is_surface_interaction());
            return (const SurfaceInteraction &)*this;
        }
    };

    // 与几何体表面的交互信息
    struct SurfaceInteraction : public Interaction {
        Vector3f dpdu {}, dpdv {};
        Normal3f dndu {}, dndv {};
        
        struct {
            Normal3f normal {};
            Vector3f dpdu {}, dpdv {};
            Normal3f dndu {}, dndv {};
        } shading {};

        int face_index = 0;

        SurfaceInteraction(
            const Point3fi &pointi, const Point2f &uv, const Vector3f &wo,
            const Vector3f &dpdu, const Vector3f &dpdv,
            const Normal3f &dndu, const Normal3f &dndv,
            Real time, bool flip_normal
        ) : Interaction(pointi, Normal3f(normalize(cross(dpdu, dpdv))), uv, wo, time), 
            dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv) {
            shading.normal = normal;
            shading.dpdu = dpdu;
            shading.dpdv = dpdv;
            shading.dndu = dndu;
            shading.dndv = dndv;
            if (flip_normal) {
                this->normal *= -1;
                shading.normal *= -1;
            }
        }
 
        void set_shading_geometry(Normal3f ns, Vector3f dpdus, Vector3f dpdvs, Normal3f dndus, Normal3f dndvs, bool orientation_is_authoritative) {
            shading.normal = ns;
            if (orientation_is_authoritative) {
                normal = face_forward(normal, shading.normal);
            } else {
                shading.normal = face_forward(shading.normal, normal);
            }
            shading.dpdu = dpdus;
            shading.dpdv = dpdvs;
            shading.dndu = dndus;
            shading.dndv = dndvs;
        }
    };

    // 与传播介质的交互信息
    struct MediumInteraction : public Interaction {
        PhaseFunction phase;

        MediumInteraction(Point3f point, Vector3f wo, Real time, Medium medium, PhaseFunction phase) : Interaction(point, wo, time, medium), phase(phase) {}
    };
}
