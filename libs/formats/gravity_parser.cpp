#include "common.h"

namespace orb {
std::string GravitationalModelParser::replaceDwithE(const std::string& s) {
    std::string result = s;
    std::replace(result.begin(), result.end(), 'D', 'E');
    std::replace(result.begin(), result.end(), 'd', 'E');
    return result;
}

GravitationalModelParser::GravitationalModelParser(const std::string& file) : file_(file) {
    parse();
}

void GravitationalModelParser::parse() {
    std::ifstream infile(file_);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open file: " + file_);
    }

    std::string line;

    // First line: Max_Degree
    if (!std::getline(infile, line))
        throw std::runtime_error("File too short");
    Max_Degree = std::stoi(line);

    // Second line: Radius
    if (!std::getline(infile, line))
        throw std::runtime_error("File too short");
    Radius = std::stod(replaceDwithE(line));

    // Third line: mu
    if (!std::getline(infile, line))
        throw std::runtime_error("File too short");
    mu = std::stod(replaceDwithE(line));

    // Parse the rest as data
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row;
        std::string token;
        while (iss >> token) {
            row.push_back(std::stod(replaceDwithE(token)));
        }
        if (!row.empty())
            data.push_back(row);
    }
}
}