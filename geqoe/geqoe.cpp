#include "geqoe.h"
#include <math.h>
#include <memory>
#include <vector>
#include "../libs/common/include/geometry.h"
#include "../libs/common/include/algo.h"
#include "perturbations.h"
#include "../libs/common/include/common.h"
#include "../libs/formats/common.h"

extern "C"{
    #include "../libs/cspice/include/SpiceUsr.h"
}

namespace orb {

    geqoe operator+(const geqoe &a, const geqoe &b)
    {
        return {a.nu + b.nu, a.p1 + b.p1, a.p2 + b.p2, a.varL + b.varL, a.q1 + b.q1, a.q2 + b.q2};
    }

    geqoe operator*(double scalar, const geqoe &a)
    {
        return {scalar * a.nu, scalar * a.p1, scalar * a.p2, scalar * a.varL, scalar * a.q1, scalar * a.q2};
    }

    GEqOE::GEqOE() = default;
    GEqOE::~GEqOE() = default;

    void GEqOE::compute(std::vector<double> tspan, geqoe s0, double dt, std::map<double, geqoe> &sprop) {
        logger.info("Starting GEqOE propagation from t = " + std::to_string(tspan[0]) + " to t = " + std::to_string(tspan[1]) + " with dt = " + std::to_string(dt));
        double t = tspan[0];
        geqoe s = s0;

        sprop[t] = s;

        while (t < tspan[1]) {
            s = runge_kutta5(s, t, dt, [this](double t, geqoe s) { return this->odes(t, s); });
            t += dt;
            sprop[t] = s;
        }
        logger.info("GEqOE propagation completed.");

        
    }

    geqoe GEqOE::odes(double t, geqoe s){
        geqoe dsdt;

        double g2 = std::pow(s.p1, 2) + std::pow(s.p2, 2);
        double a = std::pow(mu / std::pow(s.nu, 2), 1.0 / 3.0);
        double c = std::pow(std::pow(mu, 2) / s.nu, 1.0 / 3.0) * std::sqrt(1 - g2);
        double rho = a * (1 - g2);
        double alpha = 1.0 / (1.0 + std::sqrt(1.0 - g2));

        double K = solve_kepler_newton(s.p1, s.p2, s.varL, 0.0, 1e-20, 200);
        double r = a * (1 - s.p1 * std::sin(K) - s.p2 * std::cos(K));
        double r_dot = std::sqrt(mu * a) / r * (s.p2 * std::sin(K) - s.p1 * std::cos(K));

        double sL = a / r * (alpha * s.p1 * s.p2 * std::cos(K) + (1 - alpha * std::pow(s.p2, 2)) * std::sin(K) - s.p1);
        double cL = a / r * (alpha * s.p1 * s.p2 * std::sin(K) + (1 - alpha * std::pow(s.p1, 2)) * std::cos(K) - s.p2);
        double L = std::fmod(std::atan2(sL, cL), 2 * M_PI);

        double Kq = 1.0f / (1.0f + s.q1 * s.q1 + s.q2 * s.q2);

        // ex vector
        Vector3D ex;
        Vector3D ey;
        Vector3D ez;

        ex.i = Kq * (1.0f - s.q1 * s.q1 + s.q2 * s.q2);
        ex.j = Kq * (2.0f * s.q1 * s.q2);
        ex.k = Kq * (-2.0f * s.q1);

        ey.i = Kq * (2.0f * s.q1 * s.q2);
        ey.j = Kq * (1.0f + s.q1 * s.q1 - s.q2 * s.q2);
        ey.k = Kq * (2.0f * s.q2);

        ez = Vector3D::cross(ex, ey);

        Vector3D er, ef, eh;
        er = ex * std::cos(L) + ey * std::sin(L);
        ef = ey * std::cos(L) - ex * std::sin(L);
        eh = ez;

        // Position components in equinoctial frame
        double x = r * std::cos(L);
        double y = r * std::sin(L);

        double xform[6][6], xform_inv[6][6];
        sxform_c(reference_frame_.inertial.c_str(), reference_frame_.body_fixed.c_str(), t, xform);
        sxform_c(reference_frame_.body_fixed.c_str(), reference_frame_.inertial.c_str(), t, xform_inv);

        Vector3D R = x * ex + y * ey;

        Vector3D r_rot;
        r_rot.i = xform[0][0] * R.i + xform[0][1] * R.j + xform[0][2] * R.k;
        r_rot.j = xform[1][0] * R.i + xform[1][1] * R.j + xform[1][2] * R.k;
        r_rot.k = xform[2][0] * R.i + xform[2][1] * R.j + xform[2][2] * R.k;

        GravityResult gravity = gravity_model->get_gravity(r_rot, t);

        double U = gravity.potential;
        Vector3D F_rot = gravity.acceleration;

        Vector3D F = Vector3D(
            xform_inv[0][0] * F_rot.i + xform_inv[0][1] * F_rot.j + xform_inv[0][2] * F_rot.k,
            xform_inv[1][0] * F_rot.i + xform_inv[1][1] * F_rot.j + xform_inv[1][2] * F_rot.k,
            xform_inv[2][0] * F_rot.i + xform_inv[2][1] * F_rot.j + xform_inv[2][2] * F_rot.k
        );

        Vector3D omega_R = Vector3D(
            xform[3][0] * R.i + xform[3][1] * R.j + xform[3][2] * R.k,
            xform[4][0] * R.i + xform[4][1] * R.j + xform[4][2] * R.k,
            xform[5][0] * R.i + xform[5][1] * R.j + xform[5][2] * R.k
        );

        Vector3D P(0.0, 0.0, 0.0);

        double Udot = Vector3D::dot(-1.0 * F_rot, omega_R);

        auto accel_map = get_third_body_accel_map();

        for (const auto &body : third_body_)
        {
            auto it = accel_map.find(body);
            if (it != accel_map.end())
            {
                P += it->second(
                    R, t,
                    reference_frame_.inertial.c_str(),
                    gravity_cfg_.mainbody);
            }
        }

        F += P;
                 
        double h = std::sqrt(std::pow(c, 2) - 2 * std::pow(r, 2) * U);

        Vector3D V = r_dot * er + h / r * ef;

        double Pr = Vector3D::dot(P, er);
        double Pf = Vector3D::dot(P, ef);
        double Fr = Vector3D::dot(F, er);
        double Fh = Vector3D::dot(F, eh);

        float omegax = x / h * Fh;
        float omegay = y / h * Fh;
        float omegh = s.q1 * omegax - s.q2 * omegay;
        float eps_dot = Udot + r_dot * Pr + h / r * Pf;

        dsdt.nu = -3.0 * std::pow( s.nu / std::pow(mu, 2.0), 1.0 / 3.0) * eps_dot;
        dsdt.p1 = s.p2 * ((h - c) / (r * r) - omegh) + (1.0 / c) * (x / a + 2.0 * s.p2) * (2.0 * U - r * Fr) + (1.0 / (c * c)) * (y * (r + rho) + r * r * s.p1) * eps_dot;
        dsdt.p2 = s.p1 * (omegh - (h - c) / (r * r)) - (1.0 / c) * (y / a + 2.0 * s.p1) * (2.0 * U - r * Fr) + (1.0 / (c * c)) * (x * (r + rho) + r * r * s.p2) * eps_dot;
        dsdt.varL = s.nu + (h - c) / (r * r) - omegh + (1.0 / c) * (1.0 / alpha + alpha * (1.0 - r / a)) * (2.0 * U - r * Fr) + r * r_dot * alpha / (mu * c) * (r + rho) * eps_dot;
        dsdt.q1 = 0.5 * omegay * (1.0 + s.q1 * s.q1 + s.q2 * s.q2);
        dsdt.q2 = 0.5 * omegax * (1.0 + s.q1 * s.q1 + s.q2 * s.q2);

        return dsdt;
    }

    posvel GEqOE::geqoe2state(geqoe s, double t){
        posvel result;

        double a = std::pow(mu / std::pow(s.nu, 2), 1.0 / 3.0);
        double g2 = std::pow(s.p1, 2) + std::pow(s.p2, 2);
        double c = std::pow(std::pow(mu, 2) / s.nu, 1.0 / 3.0) * std::sqrt(1 - g2);
        double alpha = 1.0 / (1.0 + std::sqrt(1.0 - g2));

        double K = solve_kepler_newton(s.p1, s.p2, s.varL, s.varL, 1e-20, 200);

        double r = a * (1 - s.p1 * std::sin(K) - s.p2 * std::cos(K));
        double r_dot = std::sqrt(mu * a) / r * (s.p2 * std::sin(K) - s.p1 * std::cos(K));

        double sL = a / r * (alpha * s.p1 * s.p2 * std::cos(K) + (1 - alpha * std::pow(s.p2, 2)) * std::sin(K) - s.p1);
        double cL = a / r * (alpha * s.p1 * s.p2 * std::sin(K) + (1 - alpha * std::pow(s.p1, 2)) * std::cos(K) - s.p2);
        double L = std::atan2(sL, cL);
        if (L < 0)
        {
            L += 2 * M_PI;
        }

        double Kq = 1.0 / (1 + std::pow(s.q1, 2) + std::pow(s.q2, 2));
        Vector3D ex, ey;
        ex.i = Kq * (1 - s.q1 * s.q1 + s.q2 * s.q2);
        ex.j = Kq * (2 * s.q1 * s.q2);
        ex.k = Kq * (-2 * s.q1);

        ey.i = Kq * (2 * s.q1 * s.q2);
        ey.j = Kq * (1 + s.q1 * s.q1 - s.q2 * s.q2);
        ey.k = Kq * (2 * s.q2);

        Vector3D er = ex * std::cos(L) + ey * std::sin(L);
        Vector3D ef = -1 * ex * std::sin(L) + ey * std::cos(L);

        Vector3D R = er * r;

        double xform[6][6];
        sxform_c(reference_frame_.inertial.c_str(), reference_frame_.body_fixed.c_str(), t, xform);
        Vector3D R_rot;
        R_rot = Vector3D(
            xform[0][0] * R.i + xform[0][1] * R.j + xform[0][2] * R.k,
            xform[1][0] * R.i + xform[1][1] * R.j + xform[1][2] * R.k,
            xform[2][0] * R.i + xform[2][1] * R.j + xform[2][2] * R.k);

        GravityResult gravity = gravity_model->get_gravity(R_rot, t);
        double U = gravity.potential;

        double h = std::sqrt(std::pow(c, 2) - 2 * std::pow(r, 2) * U);
        Vector3D V = r_dot * er + h / r * ef;

        result.R = R;
        result.V = V;

        return result;
    }

    geqoe GEqOE::state2geqoe(posvel rv, double t){
        geqoe result;
        Vector3D R = rv.R;
        Vector3D V = rv.V;
        oelem oe = Kepler_elements(R, V, mu);

        double r = R.norm();
        double v = V.norm();
        double h = Vector3D::cross(R, V).norm();
        double r_dot = Vector3D::dot(R, V) / r;

        double epsk = 0.5 * std::pow(v, 2) - mu / r;

        double xform[6][6];
        sxform_c(reference_frame_.inertial.c_str(), reference_frame_.body_fixed.c_str(), t, xform);
        Vector3D R_rot;
        R_rot = Vector3D(
            xform[0][0] * R.i + xform[0][1] * R.j + xform[0][2] * R.k,
            xform[1][0] * R.i + xform[1][1] * R.j + xform[1][2] * R.k,
            xform[2][0] * R.i + xform[2][1] * R.j + xform[2][2] * R.k
        );

        GravityResult gravity = gravity_model->get_gravity(R_rot, t);
        double U = gravity.potential;

        double eps = epsk + U;
        result.nu = 1.0 / mu * std::pow((-2.0 * eps), 3.0 / 2.0);

        result.q1 = std::tan(oe.inc / 2) * std::sin(oe.raan);
        result.q2 = std::tan(oe.inc / 2) * std::cos(oe.raan);

        double Kq = 1 / (1 + std::pow(result.q1, 2) + std::pow(result.q2, 2));

        Vector3D ex, ey;
        ex.i = Kq * (1 - result.q1 * result.q1 + result.q2 * result.q2);
        ex.j = Kq * (2 * result.q1 * result.q2);
        ex.k = Kq * (-2 * result.q1);

        ey.i = Kq * (2 * result.q1 * result.q2);
        ey.j = Kq * (1 + result.q1 * result.q1 - result.q2 * result.q2);
        ey.k = Kq * (2 * result.q2);

        Vector3D er = R * (1.0 / r);

        double cL = Vector3D::dot(er, ex);
        double sL = Vector3D::dot(er, ey);
        double L = std::atan2(sL, cL);

        if (L < 0)
        {
            L += 2 * M_PI;
        }

        double Ueff = std::pow(h, 2) / 2.0 / std::pow(r, 2) + U;
        double c = std::sqrt(2.0 * std::pow(r, 2) * Ueff);
        double rho = std::pow(c, 2) / mu;
        double a = -mu / 2.0 / eps;

        result.p1 = (rho / r - 1) * std::sin(L) - c * r_dot / mu * std::cos(L);
        result.p2 = (rho / r - 1) * std::cos(L) + c * r_dot / mu * std::sin(L);

        double w = std::sqrt(mu / a);
        double sK = (mu + c * w - r * std::pow(r_dot, 2)) * std::sin(L) - r_dot * (c + w * r) * std::cos(L);
        double cK = (mu + c * w - r * std::pow(r_dot, 2)) * std::cos(L) + r_dot * (c + w * r) * std::sin(L);
        double K = std::atan2(sK, cK);

        if (K < 0)
        {
            K += 2 * M_PI;
        }

        result.varL = K + 1 / (mu + c * w) * (cK * result.p1 - sK * result.p2);

        return result;
    }

    void GEqOE::init_geqoe(std::vector<Perturbation> &perturbations, GravityConfig &gravity_cfg, ReferenceFrame &reference_frame, std::vector<ThirdBody> &third_bodies) {

        std::string gravity_model_name = "config/" + gravity_cfg.model + ".txt";

        int gravity_harmonics = 0;
        int third_body_perturbation = 0;
        for (const auto& perturbation : perturbations) {
            if (perturbation == Perturbation::GRAV_HARMONICS) {
                gravity_harmonics = 1;
            }
            else if (perturbation == Perturbation::THIRD_BODY) {
                third_body_perturbation = 1;
            }
        }

        int grav_degree = gravity_cfg.degree * gravity_harmonics;
        int grav_order = gravity_cfg.order * gravity_harmonics;

        gravity_model = std::make_unique<Gravity>(gravity_model_name, grav_degree, grav_order);

        third_body_ = third_bodies;

        if (third_body_perturbation == 0) {
            third_body_.clear();
        }

        gravity_cfg_ = gravity_cfg;
        reference_frame_ = reference_frame;
        mu = gravity_model->getMu(); // main body gravitational constant
    }

    std::unordered_map<ThirdBody, AccelFunc> GEqOE::get_third_body_accel_map()
    {
        return {
            {ThirdBody::MOON, moon_thirdbody_accel},
            {ThirdBody::SUN, sun_thirdbody_accel},
            {ThirdBody::EARTH, earth_thirdbody_accel}
        };
    }
}