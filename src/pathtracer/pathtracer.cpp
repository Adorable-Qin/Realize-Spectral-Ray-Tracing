#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"
#include <random>
#include <chrono>

#define SAMPLE_PER_COLOR 16
#define COLOR_TEMPERATURE 40000

using namespace CGL::SceneObjects;

namespace CGL {

    PathTracer::PathTracer() {
        gridSampler = new UniformGridSampler2D();
        hemisphereSampler = new UniformHemisphereSampler3D();

        tm_gamma = 2.2f;
        tm_level = 1.0f;
        tm_key = 0.18;
        tm_wht = 5.0f;
    }

    PathTracer::~PathTracer() {
        delete gridSampler;
        delete hemisphereSampler;
    }

    void PathTracer::set_frame_size(size_t width, size_t height) {
        sampleBuffer.resize(width, height);
        sampleCountBuffer.resize(width * height);
    }

    void PathTracer::clear() {
        bvh = nullptr;
        scene = nullptr;
        camera = nullptr;
        sampleBuffer.clear();
        sampleCountBuffer.clear();
        sampleBuffer.resize(0, 0);
        sampleCountBuffer.resize(0, 0);
    }

    void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                          size_t y0, size_t x1, size_t y1) {
        sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
    }

    Vector3D
    PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                    const Intersection &isect) {
        // Estimate the lighting from this intersection coming directly from a light.
        // For this function, sample uniformly in a hemisphere.

        // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
        // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

        // make a coordinate system for a hit point
        // with N aligned with the Z direction.
        Matrix3x3 o2w;
        make_coord_space(o2w, isect.n);
        Matrix3x3 w2o = o2w.T();

        // w_out points towards the source of the ray (e.g.,
        // toward the camera if this is a primary ray)
        const Vector3D hit_p = r.o + r.d * isect.t;
        const Vector3D w_out = w2o * (-r.d);

        // This is the same number of total samples as
        // estimate_direct_lighting_importance (outside of delta lights). We keep the
        // same number of samples for clarity of comparison.
        int num_samples = scene->lights.size() * ns_area_light;
        Vector3D L_out;

        // TODO (Part 3): Write your sampling loop here

        for (int i = 0; i < num_samples; i++) {
            auto w_in = hemisphereSampler->get_sample();
            auto f = isect.bsdf->f(w_out, w_in, r.wavelength) * (2 * PI);
            auto wi = o2w * w_in;
            auto nextRay = Ray(hit_p, wi);
            nextRay.min_t = EPS_F;
            nextRay.max_t = INF_F - EPS_F;
            nextRay.color = r.color;
            nextRay.wavelength = r.wavelength;

            Intersection nextIsect{};
            if (!bvh->intersect(nextRay, &nextIsect)) continue;
            L_out += f * nextIsect.bsdf->get_emission() / (double) num_samples;
        }

        return L_out;
    }

    Vector3D
    PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                    const Intersection &isect) {
        // Estimate the lighting from this intersection coming directly from a light.
        // To implement importance sampling, sample only from lights, not uniformly in
        // a hemisphere.

        // make a coordinate system for a hit point
        // with N aligned with the Z direction.
        Matrix3x3 o2w;
        make_coord_space(o2w, isect.n);
        Matrix3x3 w2o = o2w.T();

        // w_out points towards the source of the ray (e.g.,
        // toward the camera if this is a primary ray)
        const Vector3D hit_p = r.o + r.d * isect.t;
        const Vector3D w_out = w2o * (-r.d);
        Vector3D L_out;

        for (auto light: scene->lights) {
            if (!light->is_delta_light()) {
                for (int i = 0; i < ns_area_light; i++) {
                    Vector3D wi;
                    double distToLight, pdf;

                    Vector3D lightIntensity = light->sample_L(hit_p, &wi, &distToLight, &pdf);

                    if (pdf == 0) continue;

                    auto f = isect.bsdf->f(w_out, w2o * wi, r.wavelength);
                    auto nextRay = Ray(hit_p, wi);
                    nextRay.min_t = EPS_F;
                    nextRay.max_t = distToLight - EPS_F;
                    nextRay.color = r.color;
                    nextRay.wavelength = r.wavelength;

                    if (bvh->has_intersection(nextRay)) continue;
                    L_out +=
                            f * lightIntensity / pdf / (double) ns_area_light;
                }
            }
            else {
                Vector3D wi;
                double distToLight, pdf;

                Vector3D lightIntensity = light->sample_L(hit_p, &wi, &distToLight, &pdf);

                if (pdf == 0) continue;

                auto f = isect.bsdf->f(w_out, w2o * wi, r.wavelength);
                auto nextRay = Ray(hit_p, wi);
                nextRay.min_t = EPS_F;
                nextRay.max_t = distToLight - EPS_F;

                if (bvh->has_intersection(nextRay)) continue;
                L_out += f * lightIntensity / pdf;
            }
        }

        return L_out;
    }

    Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                              const Intersection &isect) {
        // TODO: Part 3, Task 2
        // Returns the light that results from no bounces of light
        // (e.g. the light emmitted from a light source)

        return isect.bsdf->get_emission();// * dot(isect.n, r.d.unit())

    }

    Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                             const Intersection &isect) {
        // TODO: Part 3, Task 3
        // Returns either the direct illumination by hemisphere or importance sampling
        // depending on `direct_hemisphere_sample`

        if (direct_hemisphere_sample)
            return estimate_direct_lighting_hemisphere(r, isect);
        else
            return estimate_direct_lighting_importance(r, isect);
    }

    Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                      const Intersection &isect) {
        Matrix3x3 o2w;
        make_coord_space(o2w, isect.n);
        Matrix3x3 w2o = o2w.T();

        Vector3D hit_p = r.o + r.d * isect.t;
        Vector3D w_out = w2o * (-r.d);

        Vector3D L_out(0, 0, 0);

        // TODO: Part 4, Task 2
        // Returns the one bounce radiance + radiance from extra bounces at this point.
        // Should be called recursively to simulate extra bounces.
        if (!isect.bsdf->is_delta())
            L_out += one_bounce_radiance(r, isect);
        // enable it if indirect only
        // if (r.depth == max_ray_depth) L_out = Vector3D(0, 0, 0);
        if (r.depth == 0) return L_out;


        Vector3D w_in;
        double pdf;

        auto f = isect.bsdf->sample_f(w_out, &w_in, &pdf, r.wavelength);
        auto wi = o2w * w_in;
        auto ray = Ray(hit_p, wi);
        ray.depth = r.depth - 1;
        ray.min_t = EPS_F;
        ray.max_t = INF_D - EPS_F;
        ray.color = r.color;
        ray.wavelength = r.wavelength;

        Intersection shadowIsect;
        if (bvh->intersect(ray, &shadowIsect) && coin_flip(0.7)) {
            auto L = at_least_one_bounce_radiance(ray, shadowIsect);
            if (isect.bsdf->is_delta())
                L += zero_bounce_radiance(ray, shadowIsect);
            auto f_l = f * L;
            L_out += f_l * abs_cos_theta(w_in) / pdf / 0.7;
        }

        return L_out;
    }

    Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
        Intersection isect;
        Vector3D L_out;

        // You will extend this in assignment 3-2.
        // If no intersection occurs, we simply return black.
        // This changes if you implement hemispherical lighting for extra credit.

        // The following line of code returns a debug color depending
        // on whether ray intersection with triangles or spheres has
        // been implemented.
        //
        // REMOVE THIS LINE when you are ready to begin Part 3.

        if (!bvh->intersect(r, &isect))
            return envLight ? envLight->sample_dir(r) : L_out;

        // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

        // TODO (Part 3): Return the direct illumination.
        L_out = zero_bounce_radiance(r, isect);
        // L_out += one_bounce_radiance(r, isect);

        // TODO (Part 4): Accumulate the "direct" and "indirect"
        // parts of global illumination into L_out rather than just direct
        L_out += at_least_one_bounce_radiance(r, isect);

        return L_out;
    }

    void PathTracer::raytrace_pixel(size_t x, size_t y) {
        // TODO (Part 1.2):
        // Make a loop that generates num_samples camera rays and traces them
        // through the scene. Return the average Vector3D.
        // You should call est_radiance_global_illumination in this function.

        // TODO (Part 5):
        // Modify your implementation to include adaptive sampling.
        // Use the command line parameters "samplesPerBatch" and "maxTolerance"

//        int num_samples = ns_aa;          // total samples to evaluate
        Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

        // My code Start

        Vector3D radiance(0.0);

        size_t num_samples = 0;

        double s1 = 0, s2 = 0, miu, sigma;

        do {
            auto rayList = std::list<Ray>();

            auto sample = origin + gridSampler->get_sample();

            // take a few samples for each color channel
            for (int i = 0; i < SAMPLE_PER_COLOR * 3; i++) {
                auto r = camera->generate_ray(sample.x / sampleBuffer.w, sample.y / sampleBuffer.h, i % 3);
                r.depth = max_ray_depth;
                rayList.push_back(r);
            }

            // combine three channels and divide by the number of channels
            auto newRadiance = Vector3D();
            for (const auto &r: rayList) {
                auto tmp = est_radiance_global_illumination(r);
                newRadiance[r.color] += tmp[r.color];
            }
            newRadiance /= (rayList.size() / 3);

            radiance = (radiance * num_samples + newRadiance) / (num_samples + 1);

            num_samples++;

            s1 += newRadiance.illum();
            s2 += newRadiance.illum() * newRadiance.illum();

            if (num_samples % samplesPerBatch == 0) {
                miu = s1 / num_samples;
                sigma = sqrt((s2 - s1 / num_samples * s1) / (num_samples - 1));

                if (1.96 * sigma < maxTolerance * miu * sqrt(num_samples))
                    break;
            }

        } while (num_samples < ns_aa);

        auto temperature = COLOR_TEMPERATURE;
        for (int color = 0; color < 3; color++) {
            auto radiance_cof = color_temperature(temperature, color);
            radiance_cof /= 255;
            radiance[color] *= radiance_cof;
        }

        sampleBuffer.update_pixel(radiance, x, y);
        sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;

        // My code End


//        sampleBuffer.update_pixel(Vector3D(0.2, 1.0, 0.8), x, y);
//        sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;

    }

    void PathTracer::autofocus(Vector2D loc) {
        Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h, 0);
        Intersection isect;

        bvh->intersect(r, &isect);

        camera->focalDistance = isect.t;
    }

    // color temperature range from 1000K to 40000K
    double PathTracer::color_temperature(int temperature, int color) {
        double radiance;
        temperature /= 100;
        switch (color)
        {
            case 0: // R
                if (temperature <= 66) {
                    radiance = 255;
                } else {
                    radiance = temperature - 60;
                    radiance = 329.698727446 * (pow(radiance, -0.1332047592));
                    if (radiance < 0) radiance = 0;
                    if (radiance > 255) radiance = 255;
                }
                break;

            case 1: // G
                if (temperature <= 66) {
                    radiance = temperature;
                    radiance = 99.4708025861 * log(radiance) - 161.1195681661;
                    if (radiance < 0) radiance = 0;
                    if (radiance > 255) radiance = 255;
                } else {
                    radiance = temperature - 60;
                    radiance = 288.1221695283 * pow(radiance, -0.0755148492);
                    if (radiance < 0) radiance = 0;
                    if (radiance > 255) radiance = 255;
                }
                break;

            case 2: // B
                if (temperature >= 66) {
                    radiance = 255;
                } else {
                    if (temperature <= 19) {
                        radiance = 0;
                    } else {
                        radiance = temperature - 10;
                        radiance = 138.5177312231 * log(radiance) - 305.0447927307;
                        if (radiance < 0) radiance = 0;
                        if (radiance > 255) radiance = 255;
                    }
                }
        }
        return radiance;
    }

} // namespace CGL
