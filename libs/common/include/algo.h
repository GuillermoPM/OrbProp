#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include "common.h"
#include "geometry.h"

namespace orb{

template<typename State, typename DerivFunc>
inline State runge_kutta5(
    const State y, double t, double dt, DerivFunc odefun
) {
    // coefficients
    // Dormand-Prince coefficients (double division!)
    const double c2 = 1.0 / 5.0;
    const double c3 = 3.0 / 10.0;
    const double c4 = 4.0 / 5.0;
    const double c5 = 8.0 / 9.0;
    const double c6 = 1.0;
    const double c7 = 1.0;

    const double a21 = 1.0 / 5.0;

    const double a31 = 3.0 / 40.0;
    const double a32 = 9.0 / 40.0;

    const double a41 = 44.0 / 45.0;
    const double a42 = -56.0 / 15.0;
    const double a43 = 32.0 / 9.0;

    const double a51 = 19372.0 / 6561.0;
    const double a52 = -25360.0 / 2187.0;
    const double a53 = 64448.0 / 6561.0;
    const double a54 = -212.0 / 729.0;

    const double a61 = 9017.0 / 3168.0;
    const double a62 = -355.0 / 33.0;
    const double a63 = 46732.0 / 5247.0;
    const double a64 = 49.0 / 176.0;
    const double a65 = -5103.0 / 18656.0;

    const double a71 = 35.0 / 384.0;
    const double a72 = 0.0;
    const double a73 = 500.0 / 1113.0;
    const double a74 = 125.0 / 192.0;
    const double a75 = -2187.0 / 6784.0;
    const double a76 = 11.0 / 84.0;

    State k1 = odefun(t, y);
    State k2 = odefun(t + c2 * dt, y + dt * (a21 * k1));
    State k3 = odefun(t + c3 * dt, y + dt * (a31 * k1 + a32 * k2));
    State k4 = odefun(t + c4 * dt, y + dt * (a41 * k1 + a42 * k2 + a43 * k3));
    State k5 = odefun(t + c5 * dt, y + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
    State k6 = odefun(t + c6 * dt, y + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));
    State k7 = odefun(t + c7 * dt, y + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6));

    return y + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6);
}


inline double solve_kepler_newton(double p1, double p2, double varL, double K0, double tol, int max_iter) {
    double K = K0;
    for (int i = 0; i < max_iter; ++i) {
        double f = K + p1 * std::cos(K) - p2 * std::sin(K) - varL;
        double df = 1.0 - p1 * std::sin(K) - p2 * std::cos(K);
        double K_new = K - f / df;
        if (std::abs(K_new - K) < tol)
            return K_new;
        K = K_new;
    }
    return K;
}


oelem Kepler_elements(Vector3D R, Vector3D V, double mu);

posvel inertial2rot(const double t, const posvel &rv);
posvel rot2inertial(const double t, const posvel &rv);
}