#include "datarecorder.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

namespace orb {
    OrbRecorder::OrbRecorder(){}

    OrbRecorder::~OrbRecorder(){}

    void OrbRecorder::record_cartesian(timeposvel &posvel_data, const std::string &filename)
    {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        file << std::setw(15) << "Time"
             << std::setw(15) << "R_x"
             << std::setw(15) << "R_y"
             << std::setw(15) << "R_z"
             << std::setw(15) << "V_x"
             << std::setw(15) << "V_y"
             << std::setw(15) << "V_z"
             << "\n";

        for (auto &[time, posvel] : posvel_data) {
            file << std::setw(15) << std::fixed << std::setprecision(6) << time
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.R.i
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.R.j
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.R.k
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.V.i
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.V.j
                 << std::setw(15) << std::fixed << std::setprecision(6) << posvel.V.k
                 << "\n"; 
        }

        file.close();
    }

    void OrbRecorder::record_geqoe(timegeqoe &geqoe_data, const std::string &filename)
    {
        std::ofstream file(filename);

        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        file << std::setw(15) << "Time"
             << std::setw(15) << "nu"
             << std::setw(15) << "p1"
             << std::setw(15) << "p2"
             << std::setw(15) << "varL"
             << std::setw(15) << "q1"
             << std::setw(15) << "q2"
             << "\n";
        file << std::string(40, '-') << std::endl;

        for (auto &[time, geqoe] : geqoe_data)
        {
            file << std::setw(15) << std::fixed << std::setprecision(6) << time
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.nu
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.p1
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.p2
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.varL
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.q1
                 << std::setw(15) << std::fixed << std::setprecision(6) << geqoe.q2
                 << "\n";
        }
        file.close();
    }

    void OrbRecorder::record_oelem(timeoelem &oelem_data, const std::string &filename) {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }
        
        file << std::setw(15) << "Time"
             << std::setw(15) << "a"
             << std::setw(15) << "e"
             << std::setw(15) << "inc"
             << std::setw(15) << "raan"
             << std::setw(15) << "argp"
             << std::setw(15) << "ta"
             << std::setw(15) << "ma"
             << "\n";
        file << std::string(40, '-') << std::endl;

        for (auto &[time, oelem] : oelem_data) {
            file << std::setw(15) << std::fixed << std::setprecision(6) << time
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.a
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.e
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.inc
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.raan
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.argp
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.ta
                 << std::setw(15) << std::fixed << std::setprecision(6) << oelem.ma
                 << "\n";
        }
        file.close();


    }
}