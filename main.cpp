#include "rgeqoe/rgeqoe.h"
#include "geqoe/geqoe.h"
#include "libs/formats/scenario_parser.h"
#include "libs/common/include/common.h"
#include "libs/common/include/datarecorder.h"
#include "libs/common/include/algo.h"
#include <fstream>
#include "libs/formats/common.h"

int main(int argc, char* argv[])
{

    if (argc < 2)
    {
        std::cerr << "Usage: ./orb <scenario.dat>\n";
        return 1;
    }


    orb::posvel rv0_test;
    orb::geqoe geqoe0, geodes;
    std::map<double, orb::geqoe> sprop;
    orb::timeposvel rvprop;

    orb::OrbRecorder recorder;

    orb::Config cfg = orb::parse_yaml_config(argv[1]);

    orb::GEqOE propagator;
    propagator.init_geqoe();

    // Convert from orbital elements to position and velocity

    orb::oelem keplerelem;
    keplerelem.a = cfg.init_cond.sma;
    keplerelem.e = cfg.init_cond.ecc;
    keplerelem.inc = cfg.init_cond.inc;
    keplerelem.raan = cfg.init_cond.raan;
    keplerelem.argp = cfg.init_cond.argp;
    keplerelem.ta = cfg.init_cond.nu;
    orb::posvel rv0 = orb::oelem2state(keplerelem, propagator.get_mu());



    std::string scenario_file = argv[1];

    // double t0 = 6.279457891829383e+08;
    // double tf = 6.362340691853830e+08;
    double t0 = cfg.time_span.start;
    double tf = cfg.time_span.end;


    orb::load_kernels("libs/cspice/kernels");


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
