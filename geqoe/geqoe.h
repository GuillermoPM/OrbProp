#pragma once
#include "../libs/common/include/geometry.h"
#include "../libs/common/include/common.h"
#include "../libs/perturbations/perturbations.h"
#include <vector>
#include <map>


namespace orb{
    geqoe operator+(const geqoe &a, const geqoe &b);
    geqoe operator*(double scalar, const geqoe &a);

    class GEqOE {

        public:
            GEqOE();
            ~GEqOE();
            Logger logger{"geqoe.log"};

            std::unique_ptr<Gravity> gravity_model; // Gravity perturbation model

            /*
            @brief Computes the propagation of the orbit using GEqOE formulation
            @param tspan Vector of time points for propagation.
            @param s0 Initial geqoe state.
            @param dt Time step for integration.
            @param sprop Output: propagated geqoe state.
            */
            void compute(std::vector<double> tspan, geqoe s0, double dt, std::map<double, geqoe> &sprop);

            /*
            @brief Computes the time derivatives of the GEqOE state.
            @param t Current time.
            @param s Current GEqOE state.
            @return Time derivatives of the GEqOE state.
            */
            geqoe odes(double t, geqoe s);

            /*
            @brief Initializes the GEqOE state from the given position and velocity vectors.
            @param rv0 Initial position and velocity vectors.
            @param t Time at which the state is initialized.
            @return The initialized GEqOE state.
            */
            posvel init_geqoe(posvel rv0, double t);

            /*
            @brief Updates the reference frame angular velocity with the given time.
            @param t Current time.
            */
            void update_av(double t);

            /*
            @brief Converts position and velocity vectors to GEqOE state.
            @param rv Position and velocity vectors.
            @param t Time at which the state is defined.
            @return The corresponding GEqOE state.
            */
            geqoe state2geqoe(posvel rv, double t);

            /*
            @brief Converts GEqOE state to position and velocity vectors.
            @param s GEqOE state.
            @param t Time at which the state is defined.
            @return The corresponding position and velocity vectors.
            */
            posvel geqoe2state(geqoe s, double t);

        private:
            double mu; // gravitational constant of the central body

    };
}