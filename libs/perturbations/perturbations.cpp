#include "perturbations.h"
#include <cmath>
#include <filesystem>
#include "../common/include/parsers.h"
#include <utility>


extern "C"
{
#include "../cspice/include/SpiceUsr.h"
}

namespace orb {


    Gravity::Gravity(const std::filesystem::path &filename, const int nHarmonics, const int mHarmonics)
        : filename(filename), nHarmonics(nHarmonics), mHarmonics(mHarmonics),
          max_degree(std::max(nHarmonics, mHarmonics)),
          zV(), zW(), zD(), zE(), zF(), zP(), zQ(), zR(),
          mu(), Radius()
    {
        initialize_coefficients();
    }
          
    Gravity::~Gravity() = default;

    GravityResult Gravity::get_gravity(const Vector3D R, const double t)
    {
        return compute_cunningham(R, nHarmonics, mHarmonics) - compute_cunningham(R, 0, 0);
        // return compute_cunningham(R, nHarmonics, mHarmonics);
    }

    
    void Gravity::initialize_coefficients() {
        GravitationalModelParser GravityModel(filename);

        // retrieve data from parser
        mu = GravityModel.getMu();
        Radius = GravityModel.getRadius();
        max_degree = GravityModel.getMaxDegree();
        data = GravityModel.getData();

        if (nHarmonics > max_degree || mHarmonics > max_degree) {
            throw std::runtime_error("Requested harmonics exceed maximum degree in the gravity model.");
        }
        

        std::vector<int> n, m;

        std::vector<double> c, s;
        for (const auto &row : data)
        {
            if (row.size() < 4)
                continue;                          
            n.push_back(static_cast<int>(row[0])); // col 0
            m.push_back(static_cast<int>(row[1])); // col 1
            c.push_back(row[2]);                   // col 2
            s.push_back(row[3]);                   // col 3
        }

        int size = max_degree + 1;
        
        zC.resize(size, std::vector<double>(size, 0.0));
        zS.resize(size, std::vector<double>(size, 0.0));

        for (size_t i = 0; i < n.size(); ++i)
        {
            int ni = n[i];
            int mi = m[i];
            if (ni >= 0 && ni <= max_degree && mi >= 0 && mi <= max_degree)
            {
                zC[ni][mi] = c[i];
                zS[ni][mi] = s[i];
            }
        }

        zV.resize(max_degree + 3, std::vector<double>(max_degree + 2, 0.0));
        zW.resize(max_degree + 3, std::vector<double>(max_degree + 2, 0.0));
        zD.resize(max_degree + 1, 0.0);
        zE.resize(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));
        zF.resize(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));
        zP.resize(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));
        zQ.resize(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));
        zR.resize(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));

        for (int n = 0; n <= max_degree; ++n) {
            if (n - 1 >= 0 && n < static_cast<int>(zV.size()) && n < static_cast<int>(zV[n - 1].size())) {
                zV[n - 1][n] = 0.0;
                zW[n - 1][n] = 0.0;
            }
            zD[n] = std::sqrt((1.0 + (n == 0 ? 1 : 0)) * (2 * n + 3) / (2 * n + 2));

            for (int m = 0; m <= n; ++m) {
                int kr_m = (m == 0) ? 1 : 0;
                int kr_m_plus1 = ((m + 1) == 0) ? 1 : 0;
                int kr_m_minus1 = ((m - 1) == 0) ? 1 : 0;

                zE[n][m] = std::sqrt((2 * n + 3.0) * (2 * n + 1) / (n + m + 1.0) / (n - m + 1.0));
                zF[n][m] = std::sqrt((2 * n + 3.0) * (n + m) * (n - m) / (2 * n - 1.0) / (n + m + 1.0) / (n - m + 1.0));
                zP[n][m] = std::sqrt((2 * n + 1.0) * (n + m + 2.0) * (n + m + 1.0) * (2 - kr_m) / (2 * n + 3.0) / (2 - kr_m_plus1)) / 2.0 / Radius;
                zQ[n][m] = std::sqrt((2 * n + 1.0) * (n - m + 2.0) * (n - m + 1.0) * (2 - kr_m) / (2 * n + 3.0) / (2 - kr_m_minus1)) / 2.0 / Radius;
                zR[n][m] = std::sqrt((2 * n + 1.0) * (n + m + 1.0) * (n - m + 1.0) / (2 * n + 3.0)) / Radius;
            }
        }
    }

    GravityResult Gravity::compute_cunningham(const Vector3D R, const unsigned int degree,  const unsigned int order) {

        double pr = R.norm();
        double rho = Radius / pr;
        double px = R.i / pr * rho;
        double py = R.j / pr * rho;
        double pz = R.k / pr * rho;
        double pr2 = rho * rho;

        const unsigned int zonals = (order > degree) ? order : degree;

        zV[0][0] = rho;
        zW[0][0] = 0.0;

        std::vector<int> nn_idx(degree + 1);
        for (int i = 0; i <= degree; ++i) {
            nn_idx[i] = i;
        }

        for (int i = 0; i <= static_cast<int>(degree); ++i)
        {
            int idx = i + 1;
            zV[idx][idx] = zD[i] * (px * zV[i][i] - py * zW[i][i]);
            zW[idx][idx] = zD[i] * (py * zV[i][i] + px * zW[i][i]);
        }

        for (int m = 0; m <= static_cast<int>(degree); ++m)
        {
            for (int n = m; n <= static_cast<int>(degree); ++n)
            {
                if (n == 0){
                    zV[n + 1][m] = zE[n][m] * pz * zV[n][m];
                    zW[n + 1][m] = zE[n][m] * pz * zW[n][m]; 
                }
                else{
                    zV[n + 1][m] = zE[n][m] * pz * zV[n][m] - zF[n][m] * pr2 * zV[n - 1][m];
                    zW[n + 1][m] = zE[n][m] * pz * zW[n][m] - zF[n][m] * pr2 * zW[n - 1][m];
                }
            }
        }

        for (int m = 0; m < 2; ++m) {
            for (int n = static_cast<int>(degree) + 1; n <= static_cast<int>(zonals); ++n)
            {
                zV[n + 1][m] = zE[n][m] * pz * zV[n][m] - zF[n][m] * pr2 * zV[n - 1][m];
                zW[n + 1][m] = zE[n][m] * pz * zW[n][m] - zF[n][m] * pr2 * zW[n - 1][m];
            }
        }

        // initialize potential and acceleration
        double potential = 0.0;
        Vector3D acceleration(0.0, 0.0, 0.0);

        std::vector<int> n_idx;
        for (int n = 1; n <= static_cast<int>(zonals); ++n) {
            n_idx.push_back(n);
        }
        for (int idx : n_idx) {
            potential += zC[idx][0] * zV[idx][0];
        }

        if (order >= 1) {
            for (int m = 1; m <= static_cast<int>(order); ++m) {
                for (int n = m; n <= static_cast<int>(zonals); ++n) {
                    potential += zC[n][m] * zV[n][m] + zS[n][m] * zW[n][m];
                }
            }
        }

        potential += -1e6 * (mu / Radius) * (potential + rho);

        acceleration = Vector3D(0.0, 0.0, 0.0);

        for (int n = 1; n <= static_cast<int>(zonals); ++n) {
            double Cnm = zC[n][0];
            acceleration.i -= Cnm * zP[n][0] * zV[n + 1][1] * 2.0;
            acceleration.j -= Cnm * zQ[n][0] * zW[n + 1][1] * 2.0;
            acceleration.k -= Cnm * zR[n][0] * zV[n + 1][0];
        }

        for (int m = 1; m <= static_cast<int>(order); ++m) {
            for (int n = m; n <= static_cast<int>(order); ++n) {
                double Cnm = zC[n][m];
                double Snm = zS[n][m];
                double Pnm = zP[n][m];
                double Qnm = zQ[n][m];
                acceleration.i += Cnm * (-Pnm * zV[n + 1][m + 1] + Qnm * zV[n + 1][m - 1])
                                + Snm * (-Pnm * zW[n + 1][m + 1] + Qnm * zW[n + 1][m - 1]);
                acceleration.j += Cnm * (-Pnm * zW[n + 1][m + 1] - Qnm * zW[n + 1][m - 1])
                                + Snm * (Pnm * zV[n + 1][m + 1] + Qnm * zV[n + 1][m - 1]);
                acceleration.k += Cnm * (-zR[n][m] * zV[n + 1][m])
                                + Snm * (-zR[n][m] * zW[n + 1][m]);
            }
        }

        acceleration = 1e3 * (mu / Radius) * (acceleration - Vector3D(px, py, pz) * rho / Radius);

        acceleration = acceleration * 0.001;
        potential *= 0.001 * 0.001;

        GravityResult result{potential, acceleration};


        return result;


    }

    Vector3D sun_thirdbody_accel(Vector3D R, double t) {
        
        Vector3D accel;

        double posvel[6];

        const char *rotating = "ITRF93";

        Vector3D R3b;

        double target = 10;
        double et = t; 
        const char *ref = "ITRF93";
        const char *abcorr = "NONE";
        int observer = 399;
        double state[6];
        double lt;

        spkez_c(target, et, ref, abcorr, observer, state, &lt);

        Vector3D Rtb(state[0], state[1], state[2]);

        Vector3D dUd = (R - Rtb) / std::pow((R - Rtb).norm(), 3);

        Vector3D dUf = Rtb / std::pow(Rtb.norm(), 3);

        double mu = 132712440018.0;

        // Compute acceleration
        accel = -mu * (dUd + dUf);

        return accel;
    }

    Vector3D moon_thirdbody_accel(Vector3D R, double t)
    {

        Vector3D accel;

        double posvel[6];

        const char *rotating = "ITRF93";

        Vector3D R3b;

        double target = 301;
        double et = t;
        const char *ref = "ITRF93";
        const char *abcorr = "NONE";
        int observer = 399;
        double state[6];
        double lt;

        spkez_c(target, et, ref, abcorr, observer, state, &lt);

        Vector3D Rtb(state[0], state[1], state[2]);

        Vector3D dUd = (R - Rtb) / std::pow((R - Rtb).norm(), 3);

        Vector3D dUf = Rtb / std::pow(Rtb.norm(), 3);

        double mu = 4.9028000661637961E+03;

        // Compute acceleration
        accel = -mu * (dUd + dUf);

        return accel;
    }

} // namespace orb