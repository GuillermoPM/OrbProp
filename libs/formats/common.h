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


// Scenario configuration parsing

enum class DataType // types of input data
{
    OELEM,
    CARTESIAN
};

enum class Perturbation // types of perturbations
{
    GRAV_HARMONICS,
    THIRD_BODY
};

enum class ThirdBody // celestial bodies for third-body perturbations
{
    MOON,
    SUN,
    EARTH
};

struct GravityConfig
{
    std::string mainbody; // main body for gravity model (e.g., Earth, Moon, Sun)
    int degree = 0; // degree of gravity model
    int order = 0; // order of gravity model
    std::string model; // gravity model file (e.g., EGM96, EGM2008)
};

struct InitCond
{
    double sma = 0.0; // semi-major axis
    double ecc = 0.0; // eccentricity
    double inc = 0.0; // inclination
    double raan = 0.0; // right ascension of ascending node
    double argp = 0.0; // argument of perigee
    double nu = 0.0; // true anomaly
};

struct TimeSpan
{
    double start = 0.0; // start time
    double end = 0.0; // end time
    double step = 0.0; // time step
};

struct ReferenceFrame // reference frame information
{
    std::string inertial; // inertial reference frame (e.g., J2000)
    std::string body_fixed; // body-fixed reference frame (e.g., IAU_EARTH, IAU_MOON, IAU_SUN)
};
}
