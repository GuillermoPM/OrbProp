#include "rgeqoe.h"
// #include "integrators.h"
#include <math.h>
#include <memory>
#include <vector>
#include "../libs/common/include/geometry.h"
#include "../libs/common/include/algo.h"
#include "../perturbations/perturbations.h"
#include "../libs/common/include/common.h"

extern "C"{
    #include "../libs/cspice/include/SpiceUsr.h"
}

namespace orb {

    rgeqoe operator+(const rgeqoe &a, const rgeqoe &b)
    {
        return {a.nu + b.nu, a.p1 + b.p1, a.p2 + b.p2, a.varL + b.varL, a.q1 + b.q1, a.q2 + b.q2};
    }

    rgeqoe operator*(double scalar, const rgeqoe &a)
    {
        return {scalar * a.nu, scalar * a.p1, scalar * a.p2, scalar * a.varL, scalar * a.q1, scalar * a.q2};
    }

    RGEqOE::RGEqOE() = default;
    RGEqOE::~RGEqOE() = default;

    void RGEqOE::compute(std::vector<double> tspan, rgeqoe s0, double dt, std::map<double, rgeqoe> &sprop) {
        logger.info("Starting RGEqOE propagation from t = " + std::to_string(tspan[0]) + " to t = " + std::to_string(tspan[1]) + " with dt = " + std::to_string(dt));
        double t = tspan[0];
        rgeqoe s = s0;

        sprop[t] = s;

        while (t < tspan[1]) {
            s = runge_kutta5(s, t, dt, [this](double t, rgeqoe s) { return this->odes(t, s); });
            t += dt;
            sprop[t] = s;
        }
        logger.info("RGEqOE propagation completed.");

        
    }

    rgeqoe RGEqOE::odes(double t,rgeqoe s){
        rgeqoe dsdt;
        
        //
        double g2 =  std::pow(s.p1,2) + std::pow(s.p2,2);
        double a = std::pow(mu/std::pow(s.nu,2), 1.0/3.0);
        double c = std::pow(std::pow(mu, 2)/s.nu, 1.0/3.0)*std::sqrt(1-g2);
        double rho = a * (1 - g2);
        double alpha = 1.0 / (1.0 + std::sqrt(1.0 - g2));

        double K = solve_kepler_newton(s.p1, s.p2, s.varL, 0.0, 1e-20, 200);
        double r = a *(1 - s.p1 * std::sin(K) - s.p2 * std::cos(K));
        double r_dot = std::sqrt(mu*a)/r*(s.p2*std::sin(K) - s.p1*std::cos(K));

        double sL = a/r*(alpha*s.p1*s.p2*std::cos(K)+(1-alpha*std::pow(s.p2,2))*std::sin(K)-s.p1);
        double cL = a/r*(alpha*s.p1*s.p2*std::sin(K)+(1-alpha*std::pow(s.p1,2))*std::cos(K)-s.p2);
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

        Vector3D R = x * ex + y * ey;

        GravityResult gravity = gravity_model->get_gravity(R, t);

        double U = gravity.potential;
        Vector3D F = gravity.acceleration;

        Vector3D P(0.0,0.0,0.0);

        P += moon_thirdbody_accel(R, t) + sun_thirdbody_accel(R,t);


        // Calculate perturbations and inertial accelerations
        Vector3D crossprod = Vector3D::cross(omega_vect, R);
        double Uc = -0.5*std::pow(crossprod.norm(),2);

        U += Uc;

        double h = std::sqrt(std::pow(c,2) - 2*std::pow(r,2)*U);

        Vector3D a_centrifugal = -1*Vector3D::cross(omega_vect, Vector3D::cross(omega_vect, R));
        Vector3D V = r_dot*er + h/r*ef;
        Vector3D a_coriolis = -2*Vector3D::cross(omega_vect, V);

        P +=  a_coriolis;
        F +=  a_centrifugal + a_coriolis;

        double Pr = Vector3D::dot(P, er);
        double Pf = Vector3D::dot(P, ef);
        double Fr = Vector3D::dot(F, er);
        double Fh = Vector3D::dot(F, eh);

        float omegax = x/h*Fh;
        float omegay = y / h *Fh;
        float omegh = s.q1 *omegax - s.q2 *omegay;
        float eps_dot = r_dot * Pr + h / r * Pf;

        dsdt.nu = -3.0 * std::pow(s.nu / mu, 2.0) * std::pow(mu, 1.0 / 3.0) * eps_dot;
        dsdt.p1 = s.p2 * ((h - c) / (r * r) - omegh) + (1.0 / c) * (x / a + 2.0 * s.p2) * (2.0 * U - r * Fr) + (1.0 / (c * c)) * (y * (r + rho) + r * r * s.p1) * eps_dot;
        dsdt.p2 = s.p1 * (omegh - (h - c) / (r * r)) - (1.0 / c) * (y / a + 2.0 * s.p1) * (2.0 * U - r * Fr) + (1.0 / (c * c)) * (x * (r + rho) + r * r * s.p2) * eps_dot;
        dsdt.varL = s.nu + (h - c) / (r * r) - omegh + (1.0 / c) * (1.0 / alpha + alpha * (1.0 - r / a)) * (2.0 * U - r * Fr) + r * r_dot * alpha / (mu * c) * (r + rho) * eps_dot;
        dsdt.q1 = 0.5 * omegay * (1.0 + s.q1 * s.q1 + s.q2 * s.q2);
        dsdt.q2 = 0.5 * omegax * (1.0 + s.q1 * s.q1 + s.q2 * s.q2);

        return dsdt;

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        return dsdt;
    }

    posvel RGEqOE::rgeqoe2state(rgeqoe s, double t)
    {
        posvel result;

        double a = std::pow(mu/std::pow(s.nu,2), 1.0/3.0);
        double g2 = std::pow(s.p1,2) + std::pow(s.p2,2);
        double c = std::pow(std::pow(mu, 2)/s.nu,1.0/3.0)*std::sqrt(1-g2);
        double alpha = 1.0 / (1.0 + std::sqrt(1.0 - g2));

        double K = solve_kepler_newton(s.p1, s.p2, s.varL, s.varL, 1e-20, 200);

        double r = a*(1-s.p1*std::sin(K)-s.p2*std::cos(K));
        double r_dot = std::sqrt(mu*a)/r*(s.p2*std::sin(K)-s.p1*std::cos(K));

        double sL = a/r*(alpha*s.p1*s.p2*std::cos(K)+(1-alpha*std::pow(s.p2,2))*std::sin(K)-s.p1);
        double cL = a/r*(alpha*s.p1*s.p2*std::sin(K)+(1-alpha*std::pow(s.p1,2))*std::cos(K)-s.p2);
        double L = std::atan2(sL, cL);
        if (L < 0){
            L += 2 * M_PI;
        }

        double Kq = 1.0/ (1 + std::pow(s.q1,2) + std::pow(s.q2,2));
        Vector3D ex, ey;
        ex.i = Kq * (1 - s.q1 * s.q1 + s.q2 * s.q2);
        ex.j = Kq * (2 * s.q1 * s.q2);
        ex.k = Kq * (-2 * s.q1);

        ey.i = Kq * (2 * s.q1 * s.q2);
        ey.j = Kq * (1 + s.q1 * s.q1 - s.q2 * s.q2);
        ey.k = Kq * (2 * s.q2);

        Vector3D er =    ex*std::cos(L) + ey*std::sin(L);
        Vector3D ef = -1*ex*std::sin(L) + ey*std::cos(L);

        Vector3D R = er * r;
        double U = 0;
        double Uc = -0.5 * std::pow(Vector3D::cross(omega_vect, R).norm(), 2);
        U += Uc;

        double h = std::sqrt(std::pow(c,2) - 2*std::pow(r,2)*U);
        Vector3D V = r_dot*er + h/r*ef;

        result.R = R;
        result.V = V;

        return result;
    }

    rgeqoe RGEqOE::state2rgeqoe(Vector3D R, Vector3D V, double t){

        rgeqoe result;

        oelem oe = Kepler_elements(R, V, mu);

        double r = R.norm();
        double v = V.norm();
        double h = Vector3D::cross(R,V).norm();
        double r_dot = Vector3D::dot(R,V)/r;

        double Uc = -0.5 * std::pow(Vector3D::cross(omega_vect, R).norm(), 2);
        double epsk = 0.5 * std::pow(v,2) - mu/r;

        // TODO-> add gravity model

        GravityResult gravity = gravity_model->get_gravity(R, t);
        double U = gravity.potential;
        U+=Uc;

        double eps = epsk + U;
        result.nu = 1.0/mu * std::pow((-2.0*eps), 3.0/2.0);

        result.q1 = std::tan(oe.inc/2)*std::sin(oe.raan);
        result.q2 = std::tan(oe.inc/2)*std::cos(oe.raan);

        double Kq = 1/ (1 + std::pow(result.q1,2) + std::pow(result.q2,2));
        
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

        if (L < 0){
            L += 2 * M_PI;
        }

        double Ueff = std::pow(h,2)/2.0/std::pow(r,2) + U;
        double c = std::sqrt(2.0*std::pow(r,2)*Ueff);
        double rho = std::pow(c,2)/mu;
        double a = -mu/2.0/eps;

        result.p1 = (rho/r-1)*std::sin(L) - c*r_dot/mu*std::cos(L);
        result.p2 = (rho/r-1)*std::cos(L) + c*r_dot/mu*std::sin(L);

        double w = std::sqrt(mu / a);
        double sK = (mu + c * w - r * std::pow(r_dot, 2)) * std::sin(L) - r_dot * (c + w * r) * std::cos(L);
        double cK = (mu + c * w - r * std::pow(r_dot, 2)) * std::cos(L) + r_dot * (c + w * r) * std::sin(L);
        double K = std::atan2(sK, cK);

        if (K < 0){
            K += 2 * M_PI;
        }

        result.varL = K + 1 / (mu + c * w) * (cK * result.p1 - sK * result.p2);

        return result;
    }


    posvel RGEqOE::init_rgeqoe(posvel rv0,double t){


        gravity_model = std::make_unique<Gravity>("config/EGM96.txt", 40, 40);

        mu = gravity_model->getMu(); // main body gravitational constant
        double xform_av[6][6], xform[6][6];
        double av[6];
        double rot[3][3];

        const char* rotating = "ITRF93";
        const char* inertial = "J2000";

        sxform_c(rotating, inertial, t, xform_av);
        sxform_c(inertial, rotating, t, xform);
        xf2rav_c(xform_av, rot, av);
        for (int i = 0; i < 3; ++i) {
            av[i] = -av[i];
        }

        // The position vector is returned in the rotating frame
        Vector3D R_rot(
            xform[0][0] * rv0.R.i + xform[0][1] * rv0.R.j + xform[0][2] * rv0.R.k +
            xform[0][3] * rv0.V.i + xform[0][4] * rv0.V.j + xform[0][5] * rv0.V.k,
            xform[1][0] * rv0.R.i + xform[1][1] * rv0.R.j + xform[1][2] * rv0.R.k +
            xform[1][3] * rv0.V.i + xform[1][4] * rv0.V.j + xform[1][5] * rv0.V.k,
            xform[2][0] * rv0.R.i + xform[2][1] * rv0.R.j + xform[2][2] * rv0.R.k +
            xform[2][3] * rv0.V.i + xform[2][4] * rv0.V.j + xform[2][5] * rv0.V.k
        );

        Vector3D V_rot(
            xform[3][0] * rv0.R.i + xform[3][1] * rv0.R.j + xform[3][2] * rv0.R.k +
            xform[3][3] * rv0.V.i + xform[3][4] * rv0.V.j + xform[3][5] * rv0.V.k,
            xform[4][0] * rv0.R.i + xform[4][1] * rv0.R.j + xform[4][2] * rv0.R.k +
            xform[4][3] * rv0.V.i + xform[4][4] * rv0.V.j + xform[4][5] * rv0.V.k,
            xform[5][0] * rv0.R.i + xform[5][1] * rv0.R.j + xform[5][2] * rv0.R.k +
            xform[5][3] * rv0.V.i + xform[5][4] * rv0.V.j + xform[5][5] * rv0.V.k
        );

        rv0.R = R_rot;
        rv0.V = V_rot;

        omega_vect = Vector3D(av[0], av[1], av[2]);

        return rv0;
    }


    void RGEqOE::update_av(double t){
        double xform[6][6];
        double av[6];
        double rot[3][3];

        const char *rotating = "ITRF93";
        const char *target = "J2000";

        sxform_c(rotating, target, t, xform);
        xf2rav_c(xform, rot, av);

        omega_vect = Vector3D(av[0], av[1], av[2]);
    }

}
