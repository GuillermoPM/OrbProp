#include "algo.h"


namespace orb {

oelem Kepler_elements(Vector3D R, Vector3D V, double mu) {

    oelem keplerelem;
    double r = R.norm();
    double v = V.norm();

    Vector3D H = Vector3D::cross(R,V);
    double h = H.norm();
    Vector3D E = -1*R/r + 1/mu * Vector3D::cross(H,V);
    double e = E.norm();

    Vector3D u1;
    if (e != 0.0) {
        u1 = E / e;
    } else {
        u1 = Vector3D(1.0, 0.0, 0.0);
    }
    Vector3D u3 = H / h;
    Vector3D u2 = Vector3D::cross(u3, u1);

    Vector3D i(1.0, 0.0, 0.0);
    Vector3D j(0.0, 1.0, 0.0);
    Vector3D k(0.0, 0.0, 1.0);

    Vector3D uln = Vector3D::cross(k, u3);

    // Inclination
    keplerelem.inc = std::acos(std::clamp(Vector3D::dot(u3, k), -1.0, 1.0));

    // Right Ascension of Ascending Node (RAAN)
    keplerelem.raan = std::atan2(Vector3D::dot(uln,j), Vector3D::dot(uln,i));
    if (keplerelem.raan < 0)
        keplerelem.raan += 2 * M_PI;

    // Argument of Periapsis (AOP)
    keplerelem.argp = std::atan2(Vector3D::dot((-1*u2), uln), Vector3D::dot(uln, u1));
    if (keplerelem.argp < 0)
        keplerelem.argp += 2 * M_PI;

    // True Anomaly (TA)
    keplerelem.ta = std::atan2(Vector3D::dot(R, u2), Vector3D::dot(R, u1));
    if (keplerelem.ta < 0)
        keplerelem.ta += 2 * M_PI;

    double u;
    if (e != 0.0) {
        u = 2.0 * std::atan(std::sqrt((1.0 - e) / (1.0 + e)) * std::tan(keplerelem.ta / 2.0));
    } else {
        u = keplerelem.ta;
    }
    keplerelem.ma = u - e * std::sin(u);

    keplerelem.a = r / (2.0 - v * v * r / mu);

    return keplerelem;
}

posvel inertial2rot(const double t, const posvel &rv)
{
    posvel result;

    double xform[6][6];

    const char *inertial = "J2000";
    const char *rotational = "ITRF93";

    sxform_c(inertial, rotational, t, xform);

    Vector3D R(
        xform[0][0] * rv.R.i + xform[0][1] * rv.R.j + xform[0][2] * rv.R.k +
        xform[0][3] * rv.V.i + xform[0][4] * rv.V.j + xform[0][5] * rv.V.k,
        xform[1][0] * rv.R.i + xform[1][1] * rv.R.j + xform[1][2] * rv.R.k +
        xform[1][3] * rv.V.i + xform[1][4] * rv.V.j + xform[1][5] * rv.V.k,
        xform[2][0] * rv.R.i + xform[2][1] * rv.R.j + xform[2][2] * rv.R.k +
        xform[2][3] * rv.V.i + xform[2][4] * rv.V.j + xform[2][5] * rv.V.k);

    Vector3D V(
        xform[3][0] * rv.R.i + xform[3][1] * rv.R.j + xform[3][2] * rv.R.k +
        xform[3][3] * rv.V.i + xform[3][4] * rv.V.j + xform[3][5] * rv.V.k,
        xform[4][0] * rv.R.i + xform[4][1] * rv.R.j + xform[4][2] * rv.R.k +
        xform[4][3] * rv.V.i + xform[4][4] * rv.V.j + xform[4][5] * rv.V.k,
        xform[5][0] * rv.R.i + xform[5][1] * rv.R.j + xform[5][2] * rv.R.k +
        xform[5][3] * rv.V.i + xform[5][4] * rv.V.j + xform[5][5] * rv.V.k);

    result.R = R;
    result.V = V;

    return result;
}

posvel rot2inertial(const double t, const posvel &rv)
{
    posvel result;

    double xform[6][6];

    const char *inertial = "J2000";
    const char *rotational = "ITRF93";

    sxform_c(rotational, inertial, t, xform);

    Vector3D R(
        xform[0][0] * rv.R.i + xform[0][1] * rv.R.j + xform[0][2] * rv.R.k +
            xform[0][3] * rv.V.i + xform[0][4] * rv.V.j + xform[0][5] * rv.V.k,
        xform[1][0] * rv.R.i + xform[1][1] * rv.R.j + xform[1][2] * rv.R.k +
            xform[1][3] * rv.V.i + xform[1][4] * rv.V.j + xform[1][5] * rv.V.k,
        xform[2][0] * rv.R.i + xform[2][1] * rv.R.j + xform[2][2] * rv.R.k +
            xform[2][3] * rv.V.i + xform[2][4] * rv.V.j + xform[2][5] * rv.V.k
        );

    Vector3D V(
        xform[3][0] * rv.R.i + xform[3][1] * rv.R.j + xform[3][2] * rv.R.k +
            xform[3][3] * rv.V.i + xform[3][4] * rv.V.j + xform[3][5] * rv.V.k,
        xform[4][0] * rv.R.i + xform[4][1] * rv.R.j + xform[4][2] * rv.R.k +
            xform[4][3] * rv.V.i + xform[4][4] * rv.V.j + xform[4][5] * rv.V.k,
        xform[5][0] * rv.R.i + xform[5][1] * rv.R.j + xform[5][2] * rv.R.k +
            xform[5][3] * rv.V.i + xform[5][4] * rv.V.j + xform[5][5] * rv.V.k
        );

    result.R = R;
    result.V = V;

    return result;
}

} // namespace orb
