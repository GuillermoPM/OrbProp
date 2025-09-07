#pragma once
#include "geometry.h"
#include <filesystem>
#include <string>
#include <iostream>

extern "C"
{
#include "../../cspice/include/SpiceUsr.h"
}
namespace orb {

    struct rgeqoe
    {
        double nu; // generalized mean anomaly
        double p1;
        double p2;
        double varL;
        double q1;
        double q2;
    };
    

    struct oelem
    {
        double a; // semi-major axis
        double e; // eccentricity
        double inc; // inclination
        double raan; // right ascension of ascending node
        double argp; // argument of periapsis
        double ta; // true anomaly
        double ma; // mean anomaly
    };

    struct posvel
    {
        Vector3D R;
        Vector3D V;
    };



    inline void load_kernels(const char* folderpath) {
        namespace fs = std::filesystem;
        if (!fs::exists(folderpath) || !fs::is_directory(folderpath))
        {
            std::cerr << "Folder not found: " << fs::absolute(folderpath) << std::endl;
            return;
        }
        for (const auto &entry : fs::directory_iterator(folderpath))
        {
            if (entry.is_regular_file())
            {
                std::string path = entry.path().string();
                std::cout << "Loading kernel: " << path << std::endl;
                furnsh_c(path.c_str());
            }
        }
    }




}
