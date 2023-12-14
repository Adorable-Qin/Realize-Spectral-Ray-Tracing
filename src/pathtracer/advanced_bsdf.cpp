#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {
    // helper function
    Vector3D BSDF::wavelength_dependent_BSDF(double wavelength, Vector3D reflectance) {
        double red_cof = 1.0 / (wavelength - 600);
        double green_cof = 1.0 / (wavelength - 550);
        double blue_cof = 1.0 / (wavelength - 450);
        double N = red_cof + green_cof + blue_cof;
        double albedo = (red_cof * reflectance[0] + green_cof * reflectance[1] + blue_cof * reflectance[2]) / N;
        return albedo;
    }

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi, double wavelength) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf, double wavelength) {

        // TODO Assignment 7: Part 1
        // Implement MirrorBSDF

        *pdf = 1.0;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Mirror BSDF")) {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Assignment 7: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        double thetah = getTheta(h.unit());
        double cos_1 = cos(thetah);
        double tan_thetah_2 = (1 - cos_1 * cos_1) / (cos_1 * cos_1);

        double up = exp(-tan_thetah_2 / (alpha * alpha));
        double down = PI * alpha * alpha * cos_1 * cos_1 * cos_1 * cos_1;

        double ans = up / down;
        return ans;
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Assignment 7: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        Vector3D F;
        double cosval = cos(getTheta(wi.unit()));

        Vector3D Rs_up = eta * eta + k * k - 2 * eta * cosval + cosval * cosval;
        Vector3D Rs_down = eta * eta + k * k + 2 * eta * cosval + cosval * cosval;
        Vector3D Rs = Rs_up / Rs_down;

        Vector3D Rp_up = (eta * eta + k * k) * cosval * cosval - 2 * eta * cosval + 1;
        Vector3D Rp_down = (eta * eta + k * k) * cosval * cosval + 2 * eta * cosval + 1;
        Vector3D Rp = Rp_up / Rp_down;

        F = (Rs + Rp) / 2;
        return F;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi, double wavelength) {
        Vector3D ans(0.0, 0.0, 0.0);
        if (wo.z >= 0) {
            Vector3D h = ((wo + wi) / 2).unit();
            Vector3D up = F(wi) * G(wo, wi) * D(h);
            double down = 4 * dot(Vector3D(0, 0, 1.0), wo) * dot(Vector3D(0, 0, 1.0), wi);
            ans = up / down;
            auto f = wavelength_dependent_BSDF(wavelength, ans);
            ans = f * wavelength / 300;
        }
        return ans;
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf, double wavelength) {
        // TODO Assignment 7: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        //*wi = cosineHemisphereSampler.get_sample(pdf);
        Vector3D ans(0., 0., 0.);
        double r1 = sampler.get_sample()[0];
        double r2 = sampler.get_sample()[1];

        double theta = atan(sqrt(-alpha * alpha * log(1 - r1)));
        double phi = 2 * PI * r2;
        Vector3D h(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

        *wi = 2. * dot(wo, h) * h - wo;
        wi->normalize();
        if (wi->z >= 0) {
            double p_theta_up = 2 * sin(theta) * exp(-tan(theta) * tan(theta) / (alpha * alpha));
            double p_theta_down = alpha * alpha * cos(theta) * cos(theta) * cos(theta);
            double p_theta = p_theta_up / p_theta_down;

            double p_phi = 1. / (2. * PI);

            double pwh = p_theta * p_phi / sin(theta);
            double pwi = pwh / (4 * dot(*wi, h));

            *pdf = pwi;
            ans = MicrofacetBSDF::f(wo, *wi, wavelength);
        }
        else {
            *pdf = 0;
        }

        return ans;
    }

    void MicrofacetBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Micofacet BSDF")) {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi, double wavelength) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf, double wavelength) {
        // TODO Assignment 7: Part 1
        // Implement RefractionBSDF

        *pdf = 1.0;
        auto t = refract(wo, wi, ior, wavelength);
        if (!t)
            return {};

        auto eta = 1.000293 / ior;
        auto fac = 589.29 / wavelength;
        eta = eta / fac;
        if (wo.z < 0)
            eta = 1.0/eta;

        return transmittance / abs_cos_theta(*wi) / pow(eta, 2.0);
    }

    void RefractionBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Refraction BSDF")) {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi, double wavelength) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf, double wavelength) {

        // TODO Assignment 7: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        *pdf = 1.0;
        auto t = refract(wo, wi, ior, wavelength);
        if (!t) {
            reflect(wo, wi);
            return reflectance / abs_cos_theta(*wi);
        }

        auto eta = 1.000293 / ior;
        auto fac = 589.29 / wavelength;
        eta = eta / fac;
        if (wo.z < 0)
            eta = 1.0/eta;

        auto R_0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
        auto R = R_0 + (1.0 - R_0) * pow(1.0 - abs_cos_theta(wo), 5.0);

        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            *pdf = 1.0 - R;
            return (1 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2.0);
        }
    }

    void GlassBSDF::render_debugger_node() {
        if (ImGui::TreeNode(this, "Refraction BSDF")) {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D *wi) {

        // TODO Assignment 7: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.

        // consider wo already in local coordinates
        // return wi in local coordinates
        *wi = Vector3D(-wo.x, -wo.y, wo.z);

    }

    bool BSDF::refract(const Vector3D wo, Vector3D *wi, double ior, double wavelength = 589) {

        // TODO Assignment 7: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.

        *wi = Vector3D();
        auto t = dot(wo, Vector3D{0, 0, 1});

        auto eta = 1.000277 / ior;      // standard atmospheric pressure at 0 degree Celsius
        auto fac = 589.29 / wavelength; // air refractive index at 589.29 nm
        eta = eta / fac;
        if (t < 0)  // existing case
            eta = 1.0 / eta;

        auto internal_term = 1.0 - eta * eta * (1.0 - pow(wo.z, 2.0));
        if (internal_term < 0)
            return false;

        wi->x = -eta * wo.x;
        wi->y = -eta * wo.y;
        wi->z = sqrt(internal_term);

        if (t > 0)
            wi->z = -wi->z;

        return true;
    }

} // namespace CGL
