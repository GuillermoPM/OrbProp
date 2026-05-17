#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include "../common/include/geometry.h"

namespace orb {
class GravitationalModelParser
{
public:
    GravitationalModelParser(const std::string &filename);

    bool parse();

    int getMaxDegree() const;
    double getRadius() const;
    double getMu() const;

    // Returns the parsed data as a 2D vector of doubles
    const std::vector<std::vector<double>> &getData() const;

private:
    std::string file;
    int Max_Degree;
    double Radius;
    double mu;
    std::vector<std::vector<double>> data;
};

struct OrbitalInput {
    double t0; // initial time
    double tf; // final time
    Vector3D R0; // initial position
    Vector3D V0; // initial velocity
};

class ConfigParser {
    public:
        using SectionMap = std::unordered_map<std::string, double>;
        using ConfigMap  = std::unordered_map<std::string, SectionMap>;

        static ConfigMap parse(const std::string& text);
};
}
