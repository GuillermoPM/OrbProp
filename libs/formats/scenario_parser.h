#pragma once

#include <string>
#include <vector>
#include "common.h"

namespace YAML
{
    class Node;
}

namespace orb
{

    struct Config
    {
        DataType type;

        std::vector<Perturbation> perturbations;
        std::vector<ThirdBody> third_bodies;

        GravityConfig gravity;
        InitCond init_cond;
        TimeSpan time_span;
        ReferenceFrame reference_frame;
    };

    Config parse_yaml_config(const std::string &filename);

} // namespace orb