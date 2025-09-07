#include "rgeqoe/rgeqoe.h"
#include "libs/common/include/common.h"
#include "libs/common/include/algo.h"
#include <fstream>
#include "libs/formats/common.h"

int main()
{
    orb::posvel rv0, rv0_test;
    orb::rgeqoe rgeqoe0, rgeodes;
    std::map<double, orb::rgeqoe> sprop;
    std::vector<orb::posvel> rvprop;

    rv0.R = orb::Vector3D(3.063755995412349e+03, 1.174417162099421e+04, 1.780952065585180e+03);
    rv0.V = orb::Vector3D(1.701021806282415, -1.275844584985527, 5.290694112967110);
    double t0 = 6.279457891829383e+08;
    double tf = 6.362340691853830e+08;


    orb::load_kernels("libs/cspice/kernels");
    orb::RGEqOE propagator;
    orb::posvel rv0_ecef = propagator.init_rgeqoe(rv0, t0);

    rgeqoe0 = propagator.state2rgeqoe(rv0_ecef.R, rv0_ecef.V, t0);
    rv0_test = propagator.rgeqoe2state(rgeqoe0, t0);

    propagator.compute({t0, tf}, rgeqoe0, 50, sprop);

    for (const auto& s : sprop) {
        orb::posvel rv_rot, rv;
        double t = s.first;
        rgeodes = s.second;

        rv_rot = propagator.rgeqoe2state(rgeodes, t);
        rv = orb::rot2inertial(t, rv_rot);

        rvprop.push_back(rv);
    }

    std::ofstream outfile("rvprop_output.csv");
    outfile << "R_i,R_j,R_k,V_i,V_j,V_k\n";
    for (const auto& rv : rvprop) {
        outfile << rv.R.i << "," << rv.R.j << "," << rv.R.k << ","
                << rv.V.i << "," << rv.V.j << "," << rv.V.k << "\n";
    }
    outfile.close();


    return 0;
}
