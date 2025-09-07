#include "parsers.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace orb{
    ConfigParser::ConfigParser(const std::string& filename) {
    }

    bool ConfigParser::parse() {
        // Parsing implementation
        return true;
    }

    double ConfigParser::getValue(const std::string& key) const {
        // Get value implementation
        return 0.0;
    }

    GravitationalModelParser::GravitationalModelParser(const std::string& filename)
        : file(filename), Max_Degree(0), Radius(0.0), mu(0.0)
    {
        parse();
    }

    int GravitationalModelParser::getMaxDegree() const {
        return Max_Degree;
    }

    double GravitationalModelParser::getRadius() const {
        return Radius;
    }

    double GravitationalModelParser::getMu() const {
        return mu;
    }

    const std::vector<std::vector<double>>& GravitationalModelParser::getData() const {
        return data;
    }

    bool GravitationalModelParser::parse() {
        std::ifstream f(file);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open file: " + file);
        }

        auto clean_number = [](std::string s)
        {
            s.erase(std::remove_if(s.begin(), s.end(),
                                   [](unsigned char c)
                                   { return std::isspace(c); }),
                    s.end());
            std::replace(s.begin(), s.end(), 'D', 'E');
            return s;
        };
        std::string line;

        // Max degree
        std::getline(f, line);
        Max_Degree = std::stoi(clean_number(line));

        // Radius
        std::getline(f, line);
        Radius = std::stod(clean_number(line));

        // mu
        std::getline(f, line);
        mu = std::stod(clean_number(line));

        // Remaining data
        data.clear();
        while (std::getline(f, line))
        {
            std::istringstream iss(line);
            std::vector<double> row;
            std::string value;
            while (iss >> value)
            {
                row.push_back(std::stod(clean_number(value)));
            }
            if (!row.empty())
            {
                data.push_back(row);
            }
        }

        return true;
    }



}