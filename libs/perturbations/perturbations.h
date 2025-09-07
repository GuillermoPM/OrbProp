#pragma once
#include <vector>
#include <filesystem>
#include "../common/include/geometry.h"
#include <utility>

namespace orb {
    struct GravityResult {
        double potential;      // gravitational potential
        Vector3D acceleration; // gravitational acceleration
    };

    inline GravityResult operator-(const GravityResult& lhs, const GravityResult& rhs) {
        return GravityResult{
            lhs.potential - rhs.potential,
            lhs.acceleration - rhs.acceleration
        };
    }

    class Gravity {
        public:


            Gravity(const std::filesystem::path &filename, const int nHarmonics, const int mHarmonics);
            ~Gravity();

            GravityResult get_gravity(const Vector3D R, const double t);
            double getMu() const { return mu; }
            double getRadius() const { return Radius; }
            int getMaxDegree() const { return max_degree; }

            Gravity operator-(const Gravity &other) const;

        private:
            // function to compute potential and acceleration

            std::vector<std::vector<double>> zV;
            std::vector<std::vector<double>> zW;
            std::vector<double> zD;
            std::vector<std::vector<double>> zE, zF, zP, zQ, zR;
            std::vector<std::vector<double>> zC, zS;

            GravityResult compute_cunningham(const Vector3D R, const unsigned int degree, const unsigned int order);
            void initialize_coefficients();

            // filename data
            std::filesystem::path filename;
            int nHarmonics;
            int mHarmonics;
            double mu;
            double Radius;
            unsigned int max_degree;
            std::vector<std::vector<double>> data;


    };
    Vector3D thirdbody_accel(Vector3D R, double t);
} 