#pragma once

#include <MatchOfflineRenderer/camera/camera_base.hpp>
#include <MatchOfflineRenderer/sample/sample_unit_disk.hpp>
#include <MatchOfflineRenderer/math/sphere.hpp>

namespace MatchOfflineRenderer::camera {
    // 投影相机
    struct ProjectiveCamera : public Camera {
        math::Transform camera_to_screen {}, camera_to_raster {}, screen_to_raster {};
        math::Real lens_radius {};
        math::Real focal_distance {};

        ProjectiveCamera(const CameraCreateInfo &info, const math::Transform &camera_to_screen, const math::Bounds2f &screen_window, math::Real lens_radius, math::Real focal_distance) noexcept : Camera(info), camera_to_screen(camera_to_screen), lens_radius(lens_radius), focal_distance(focal_distance) {
            math::Transform screen_to_ndc = math::Transform::generate_scale(
                1 / (screen_window.max.x - screen_window.min.x),
                1 / (screen_window.max.y - screen_window.min.y),
                1
            ) * math::Transform::generate_translate({ -screen_window.min.x, -screen_window.min.y, 0 });
            math::Transform ndc_to_raster = math::Transform::generate_scale(
                film->resolution.x,
                -film->resolution.y,
                1
            );
            screen_to_raster = ndc_to_raster * screen_to_ndc;
            camera_to_raster = screen_to_raster * camera_to_screen;
        }
    };

    // 正交投影相机
    struct OrthographicCamera : public ProjectiveCamera {
        math::Vector3f origin_dx {}, origin_dy {};
        
        OrthographicCamera(const CameraCreateInfo &info, const math::Bounds2f &screen_window, math::Real lens_radius, math::Real focal_distance) noexcept : ProjectiveCamera(info, math::Transform::generate_orthographic(0, 1), screen_window, lens_radius, focal_distance) {
            origin_dx = camera_to_raster.apply_inverse(math::Vector3f { 1, 0, 0 });
            origin_dy = camera_to_raster.apply_inverse(math::Vector3f { 0, 1, 0 });
        }
        
        std::optional<CameraRay> generate_ray(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept override {
            math::Point3f raster_film_point { sample.film_point.x, sample.film_point.y, 0 };
            math::Point3f camera_film_point = camera_to_raster.apply_inverse(raster_film_point);
            math::Ray ray { camera_film_point, { 0, 0, 1 }, { sample_time(sample.time_sample_t), &medium } };

            if (lens_radius > 0) {
                auto lens_point = sample::sample_uniform_disk_concentric(sample.lens_sample_t) * lens_radius;
                math::Real ft = focal_distance / ray.direction.z;
                auto focus_point = ray(ft);
                ray.origin = { lens_point.x, lens_point.y, 0 };
                ray.direction = math::normalize(focus_point - ray.origin);
            }

            return CameraRay{ camera_transform.camera_to_render(ray) };
        }

        // std::optional<CameraRayDifferential> generate_ray_differential(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept override {
        //     math::Point3f raster_film_point { sample.film_point.x, sample.film_point.y, 0 };
        //     math::Point3f camera_film_point = camera_to_raster.apply_inverse(raster_film_point);
        //     math::RayDifferential ray_differential { camera_film_point, { 0, 0, 1 }, { sample_time(sample.time_sample_t), &medium } };

        //     if (lens_radius > 0) {
        //         auto lens_point = sample::sample_uniform_disk_concentric(sample.lens_sample_t) * lens_radius;
        //         math::Real ft = focal_distance / ray_differential.direction.z;
        //         auto focus_point = ray_differential(ft);
        //         ray_differential.origin = { lens_point.x, lens_point.y, 0 };
        //         ray_differential.direction = math::normalize(focus_point - ray_differential.origin);

        //         ray_differential.origin_x = ray_differential.origin;
        //         ray_differential.origin_y = ray_differential.origin;

        //         auto focus_point_x = focus_point + origin_dx;
        //         ray_differential.direction_x = math::normalize(focus_point_x - ray_differential.origin);
        //         auto focus_point_y = focus_point + origin_dy;
        //         ray_differential.direction_y = math::normalize(focus_point_y - ray_differential.origin);
        //     } else {
        //         ray_differential.origin_x = ray_differential.origin + origin_dx;
        //         ray_differential.origin_y = ray_differential.origin + origin_dy;
        //         ray_differential.direction_x = ray_differential.direction;
        //         ray_differential.direction_y = ray_differential.direction;
        //     }

        //     ray_differential.has_differentials = true;
        //     return CameraRayDifferential { camera_transform.camera_to_render(ray_differential) };
        // }
    };

    // 透视投影相机
    struct PerspectiveCamera : public ProjectiveCamera {
        math::Vector3f origin_dx {}, origin_dy {};

        PerspectiveCamera(const CameraCreateInfo &info, const math::Bounds2f &screen_window, math::Real fov, math::Real lens_radius, math::Real focal_distance) noexcept : ProjectiveCamera(info, math::Transform::generate_perspective(fov, 1e-2, 1000.0), screen_window, lens_radius, focal_distance) {
            origin_dx = camera_to_raster.apply_inverse(math::Vector3f { 1, 0, 0 }) - camera_to_raster.apply_inverse(math::Vector3f { 0, 0, 0 });
            origin_dy = camera_to_raster.apply_inverse(math::Vector3f { 0, 1, 0 }) - camera_to_raster.apply_inverse(math::Vector3f { 0, 0, 0 });

            MCH_ERROR("Not Implement")
            // Point2f radius = Point2f(film.GetFilter().Radius());
            // Point3f pCorner(-radius.x, -radius.y, 0.f);
            // Vector3f wCornerCamera = Normalize(Vector3f(cameraFromRaster(pCorner)));
            // cosTotalWidth = wCornerCamera.z;

            // FindMinimumDifferentials(this);
        }

        std::optional<CameraRay> generate_ray(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept override {
            math::Point3f raster_film_point { sample.film_point.x, sample.film_point.y, 0 };
            math::Point3f camera_film_point = camera_to_raster.apply_inverse(raster_film_point);
            math::Ray ray { { 0, 0, 0 }, math::normalize(math::Vector3f { camera_film_point }), { sample_time(sample.time_sample_t), &medium } };

            if (lens_radius > 0) {
                auto lens_point = sample::sample_uniform_disk_concentric(sample.lens_sample_t) * lens_radius;
                math::Real ft = focal_distance / ray.direction.z;
                auto focus_point = ray(ft);
                ray.origin = { lens_point.x, lens_point.y, 0 };
                ray.direction = math::normalize(focus_point - ray.origin);
            }

            return CameraRay{ camera_transform.camera_to_render(ray) };
        }

        // std::optional<CameraRayDifferential> generate_ray_differential(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept override {
        //     math::Point3f raster_film_point { sample.film_point.x, sample.film_point.y, 0 };
        //     math::Point3f camera_film_point = camera_to_raster.apply_inverse(raster_film_point);
        //     math::RayDifferential ray_differential { { 0, 0, 0 }, math::normalize(math::Vector3f { camera_film_point }), { sample_time(sample.time_sample_t), &medium } };

        //     if (lens_radius > 0) {
        //         auto lens_point = sample::sample_uniform_disk_concentric(sample.lens_sample_t) * lens_radius;
        //         math::Real ft = focal_distance / ray_differential.direction.z;
        //         auto focus_point = ray_differential(ft);
        //         ray_differential.origin = { lens_point.x, lens_point.y, 0 };
        //         ray_differential.direction = math::normalize(focus_point - ray_differential.origin);

        //         ray_differential.origin_x = ray_differential.origin;
        //         ray_differential.origin_y = ray_differential.origin;

        //         auto focus_point_x = focus_point + origin_dx;
        //         ray_differential.direction_x = math::normalize(focus_point_x - ray_differential.origin);
        //         auto focus_point_y = focus_point + origin_dy;
        //         ray_differential.direction_y = math::normalize(focus_point_y - ray_differential.origin);
        //     } else {
        //         ray_differential.origin_x = ray_differential.origin;
        //         ray_differential.origin_y = ray_differential.origin;
        //         ray_differential.direction_x = normalize(math::Vector3f (camera_film_point) + origin_dx);
        //         ray_differential.direction_y = normalize(math::Vector3f (camera_film_point) + origin_dy);
        //     }

        //     ray_differential.has_differentials = true;
        //     return CameraRayDifferential { camera_transform.camera_to_render(ray_differential) };
        // }
    };

    // 球相机
    struct SphericalCamera : public Camera {
        enum class Mapping {
            eEquiRectangular,
            eEqualArea,
        } mapping;

        SphericalCamera(const CameraCreateInfo &info, Mapping mapping) noexcept : Camera(info), mapping(mapping) {
            MCH_ERROR("Not Implement")
            // FindMinimumDifferentials(this);
        }
        
        std::optional<CameraRay> generate_ray(const CameraSample &sample, spectra::SampledWavelengths &wavelengths) const noexcept override {
            math::Point2f uv = { sample.film_point.x / film->resolution.x, sample.film_point.y / film->resolution.y };
            math::Vector3f direction {};

            switch (mapping) {
            case Mapping::eEquiRectangular: {
                math::Real theta = math::pi * uv.y, phi = 2 * math::pi * uv.x;
                direction = math::spherical_direction(std::sin(theta), std::cos(theta), phi);
                break;
            }
            case Mapping::eEqualArea: {
                direction = math::equal_area_square_to_sphere(math::wrap_equal_area_square(uv));
                break;
            }
            }
            std::swap(direction.y, direction.z);
            math::Ray ray { { 0, 0, 0 }, direction, { sample_time(sample.time_sample_t), &medium } };
            return CameraRay{ camera_transform.camera_to_render(ray) };
        }
    };
}
