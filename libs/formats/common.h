#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "../common/include/geometry.h"

namespace orb {
class GravitationalModelParser {
public:
    int Max_Degree;
    double Radius;
    double mu;
    std::vector<std::vector<double>> data;

    explicit GravitationalModelParser(const std::string& file);

private:
    std::string file_;

    static std::string replaceDwithE(const std::string& s);
    void parse();
};

struct OrbitalInput {
    double t0; // initial time
    double tf; // final time
    Vector3D R0; // initial position
    Vector3D V0; // initial velocity
};

class OrbParser {
    public:
        // constructor
        explicit OrbParser(const std::string& file);

        void parse();

        OrbitalInput getData() const;

    private:
        std::string file_;
        OrbitalInput data_;

        void parseLine(const std::string& line);
};


}