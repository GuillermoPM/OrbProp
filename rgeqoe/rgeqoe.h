#pragma once
#include "../libs/common/include/geometry.h"
#include "../libs/common/include/common.h"
#include "../libs/perturbations/perturbations.h"
#include <vector>
#include <map>

namespace orb {

    rgeqoe operator+(const rgeqoe &a, const rgeqoe &b);
    rgeqoe operator*(double scalar, const rgeqoe &a);

    class RGEqOE {

        public:
            RGEqOE();
            ~RGEqOE();
            Logger logger{"rgeqoe.log"};

            std::unique_ptr<Gravity> gravity_model; // Gravity perturbation model

            /*
            @brief Computes the propagation of the orbit using RGEqOE formulation
            @param tspan Vector of time points for propagation.
            @param s0 Initial rgeqoe state.
            @param dt Time step for integration.
            @param sprop Output: propagated rgeqoe state.
            */
            void compute(std::vector<double> tspan, rgeqoe s0, double dt, std::map<double, rgeqoe> &sprop);

            /*
            @brief Computes the time derivatives of the RGEqOE state.
            @param t Current time.
            @param s Current RGEqOE state.
            @return Time derivatives of the RGEqOE state.
            */
            rgeqoe odes(double t, rgeqoe s);

            /*
            @brief Initializes the RGEqOE state from the given position and velocity vectors.
            @param rv0 Initial position and velocity vectors.
            @param t Time at which the state is initialized.
            @return The initialized RGEqOE state.
            */
            posvel init_rgeqoe(posvel rv0, double t);

            /*
            @brief Updates the reference frame angular velocity with the given time.
            @param t Current time.
            */
            void update_av(double t);

            /*
            @brief Converts position and velocity vectors to RGEqOE state.
            @param R Position vector.
            @param V Velocity vector.
            @param t Time at which the state is defined.
            @return The corresponding RGEqOE state.
            */
            rgeqoe state2rgeqoe(Vector3D R, Vector3D V, double t);

            /*
            @brief Converts RGEqOE state to position and velocity vectors.
            @param s RGEqOE state.
            @param t Time at which the state is defined.
            @return The corresponding position and velocity vectors.
            */
            posvel rgeqoe2state(rgeqoe s, double t);

        private:
            double mu;              // gravitational parameter
            Vector3D omega_vect;    // reference frame angular velocity

    };



}
