#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

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
}