#include "scenario_parser.h"
#include "common.h"

#include <yaml-cpp/yaml.h>
#include <stdexcept>

namespace orb
{

    static DataType parse_data_type(const std::string &s)
    {
        if (s == "oelem")
            return DataType::OELEM;
        if (s == "cartesian")
            return DataType::CARTESIAN;

        throw std::runtime_error("Unknown TYPE: " + s);
    }

    static Perturbation parse_perturbation(const std::string &s)
    {
        if (s == "grav_harmonics")
            return Perturbation::GRAV_HARMONICS;
        if (s == "third_body")
            return Perturbation::THIRD_BODY;

        throw std::runtime_error("Unknown perturbation: " + s);
    }

    static ThirdBody parse_third_body(const std::string &s)
    {
        if (s == "Moon")
            return ThirdBody::MOON;
        if (s == "Sun")
            return ThirdBody::SUN;
        if (s == "Earth")
            return ThirdBody::EARTH;

        throw std::runtime_error("Unknown third body: " + s);
    }

    Config parse_yaml_config(const std::string &filename)
    {
        YAML::Node root = YAML::LoadFile(filename);

        Config cfg;

        // TYPE
        cfg.type = parse_data_type(root["type"].as<std::string>());

        // PERTURBATIONS
        if (root["perturbations"])
        {
            for (const auto &p : root["perturbations"])
                cfg.perturbations.push_back(parse_perturbation(p.as<std::string>()));
        }

        // GRAVITY
        if (root["gravity"])
        {
            auto g = root["gravity"];
            cfg.gravity.mainbody = g["mainbody"].as<std::string>();
            cfg.gravity.degree = g["degree"].as<int>();
            cfg.gravity.order = g["order"].as<int>();
            cfg.gravity.model = g["model"].as<std::string>();
        }

        // THIRD BODIES
        if (root["third_bodies"])
        {
            for (const auto &b : root["third_bodies"])
                cfg.third_bodies.push_back(parse_third_body(b.as<std::string>()));
        }

        // INIT COND
        if (root["initial_conditions"])
        {
            auto ic = root["initial_conditions"];

            cfg.init_cond.sma = ic["sma"].as<double>();
            cfg.init_cond.ecc = ic["ecc"].as<double>();
            cfg.init_cond.inc = ic["inc"].as<double>();
            cfg.init_cond.raan = ic["raan"].as<double>();
            cfg.init_cond.argp = ic["argp"].as<double>();
            cfg.init_cond.nu = ic["nu"].as<double>();
        }

        if (root["reference_frame"])
        {
            auto rf = root["reference_frame"];
            cfg.reference_frame.inertial = rf["inertial"].as<std::string>();
            cfg.reference_frame.body_fixed = rf["body_fixed"].as<std::string>();
        }

        // TIME SPAN
        if (root["time_span"])
        {
            auto ts = root["time_span"];

            cfg.time_span.start = ts["start"].as<double>();
            cfg.time_span.end = ts["end"].as<double>();
            cfg.time_span.step = ts["step"].as<double>();
        }

        // sanity checks
        if (cfg.gravity.order > cfg.gravity.degree)
            throw std::runtime_error("Gravity order > degree");

        if (cfg.time_span.step <= 0.0)
            throw std::runtime_error("Invalid time step");

        return cfg;
    }

} // namespace orb