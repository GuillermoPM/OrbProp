#include "rgeqoe/rgeqoe.h"
#include "geqoe/geqoe.h"
#include "libs/common/include/common.h"
#include "libs/common/include/datarecorder.h"
#include "libs/common/include/algo.h"
#include <fstream>
#include "libs/formats/common.h"

int main(int argc, char* argv[])
{
    orb::posvel rv0, rv0_test;
    orb::geqoe geqoe0, geodes;
    std::map<double, orb::geqoe> sprop;
    orb::timeposvel rvprop;

    orb::OrbRecorder recorder;

    if (argc < 2)
    {
        std::cerr << "Usage: ./orb <scenario.dat>\n";
        return 1;
    }

    std::string scenario_file = argv[1];


    rv0.R = orb::Vector3D(3.063755995412349e+03, 1.174417162099421e+04, 1.780952065585180e+03);
    rv0.V = orb::Vector3D(1.701021806282415, -1.275844584985527, 5.290694112967110);
    double t0 = 6.279457891829383e+08;
    double tf = 6.362340691853830e+08;


    orb::load_kernels("libs/cspice/kernels");
    orb::GEqOE propagator;
    propagator.init_geqoe();

    geqoe0 = propagator.state2geqoe(rv0, t0);
    rv0_test = propagator.geqoe2state(geqoe0, t0);

    propagator.compute({t0, tf}, geqoe0, 50, sprop);

    recorder.record_geqoe(sprop, "geqoe_output.txt");
    for (const auto& s : sprop) {
        orb::posvel rv_rot, rv;
        double t = s.first;
        geodes = s.second;

        rv = propagator.geqoe2state(geodes, t);

        rvprop[t] = rv;
    }


    recorder.record_cartesian(rvprop, "rvprop_output.txt");
    // recorder.record_oelem(oelem_data, "oelem_output.txt");




    return 0;
}
