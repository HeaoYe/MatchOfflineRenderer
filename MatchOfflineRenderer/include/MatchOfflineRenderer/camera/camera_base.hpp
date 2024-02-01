#pragma once

#include <MatchOfflineRenderer/camera/film.hpp>
#include <MatchOfflineRenderer/math/ray.hpp>
#include <MatchOfflineRenderer/math/transform.hpp>
#include <MatchOfflineRenderer/math/interaction.hpp>
#include <MatchOfflineRenderer/spectra/specified_spectrum.hpp>
#include <optional>

namespace MatchOfflineRenderer::camera {
    // 相机采样点
    struct CameraSample {
        math::Point2f film_point {};
        math::Point2f lens_sample_t {};
        math::Real time_sample_t {};
        math::Real filter_weight {};
    };

    // 相机光线
    struct CameraRay {
        math::Ray ray {};
        spectra::SampledSpectrum weight { 1 };
    };

    // struct CameraRayDifferential {
    //     math::RayDifferential ray {};
    //     spectra::SampledSpectrum weight { 1 };
    // };

    // 相机变换信息
    struct CameraTransform {
        math::Transform camera_to_render {};
        math::Transform render_to_world {};

        CameraTransform() noexcept = default;

        CameraTransform(const math::Transform &camera_to_world) noexcept {
            auto camera_pos = camera_to_world(math::Point3f { 0, 0, 0 });
            render_to_world = math::Transform::generate_translate(math::Vector3f(camera_pos));
            auto world_to_render = math::inverse(render_to_world);
            camera_to_render = world_to_render * camera_to_world;
        }
    };

    struct CameraCreateInfo {
        const CameraTransform &camera_transform;
        math::Real shutter_open;
        math::Real shutter_close;
        Film *film;
        const math::Medium &medium;
    };

    // 相机
    struct Camera {
        CameraTransform camera_transform {};
        math::Real shutter_open {}, shutter_close {};
        Film *film;
        math::Medium medium {};
        ImageMetadata *metadata { nullptr };

        Camera(const CameraCreateInfo &info) noexcept : camera_transform(info.camera_transform), shutter_open(info.shutter_open), shutter_close(info.shutter_close), film(info.film), medium(info.medium) {}
        
        virtual std::optional<CameraRay> generate_ray(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept = 0;
        virtual ~Camera() noexcept = default;

        virtual void setup_metadata(ImageMetadata *metadata) noexcept {
        }

        // virtual std::optional<CameraRayDifferential> generate_ray_differential(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept {
        //     auto camera_ray = generate_ray(sample, wavelengths);
        //     if (!camera_ray.has_value()) {
        //         return {};
        //     }
        //     math::RayDifferential ray_differential { camera_ray->ray };

        //     std::optional<CameraRay> camera_ray_x;
        //     for (math::Real eps : { .05f, -.05f }) {
        //         CameraSample sample_shift = sample;
        //         sample_shift.film_point.x += eps;
        //         camera_ray_x = generate_ray(sample_shift, wavelengths);
        //         if (camera_ray_x.has_value()) {
        //             ray_differential.origin_x = ray_differential.origin + (camera_ray_x->ray.origin - ray_differential.origin) / eps;
        //             ray_differential.direction_x = ray_differential.direction + (camera_ray_x->ray.direction - ray_differential.direction) / eps;
        //             break;
        //         }
        //     }
        //     std::optional<CameraRay> camera_ray_y;
        //     for (math::Real eps : { .05f, -.05f }) {
        //         CameraSample sample_shift = sample;
        //         sample_shift.film_point.y += eps;
        //         camera_ray_y = generate_ray(sample_shift, wavelengths);
        //         if (camera_ray_y.has_value()) {
        //             ray_differential.origin_y = ray_differential.origin + (camera_ray_y->ray.origin - ray_differential.origin) / eps;
        //             ray_differential.direction_y = ray_differential.direction + (camera_ray_y->ray.direction - ray_differential.direction) / eps;
        //             break;
        //         }
        //     }
        //     ray_differential.has_differentials = camera_ray_x.has_value() && camera_ray_y.has_value();
        //     return CameraRayDifferential { ray_differential, camera_ray->weight };
        // }

        virtual math::Real sample_time(math::Real t) const noexcept {
            return math::lerp(shutter_open, shutter_close, t);
        }
    };
}
