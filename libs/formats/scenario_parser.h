#pragma once

#include <string>
#include <vector>

namespace YAML
{
    class Node;
}

namespace orb
{

    enum class DataType
    {
        OELEM,
        CARTESIAN
    };

    enum class Perturbation
    {
        GRAV_HARMONICS,
        THIRD_BODY
    };

    enum class ThirdBody
    {
        MOON,
        SUN
    };

    struct GravityConfig
    {
        int degree = 0;
        int order = 0;
        std::string model;
    };

    struct InitCond
    {
        double sma = 0.0;
        double ecc = 0.0;
        double inc = 0.0;
        double raan = 0.0;
        double argp = 0.0;
        double nu = 0.0;
    };

    struct TimeSpan
    {
        double start = 0.0;
        double end = 0.0;
        double step = 0.0;
    };

    struct Config
    {
        DataType type;

        std::vector<Perturbation> perturbations;
        std::vector<ThirdBody> third_bodies;

        GravityConfig gravity;
        InitCond init_cond;
        TimeSpan time_span;
    };

    Config parse_yaml_config(const std::string &filename);

} // namespace orb