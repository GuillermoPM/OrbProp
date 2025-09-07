#pragma once
#include <string>
#include <map>
#include <vector>

namespace orb{

    class ConfigParser {
    public:
        ConfigParser(const std::string& filename);
        bool parse();
        double getValue(const std::string& key) const;

        using cfg_ = std::map<std::string, std::map<std::string, std::string>>;
    };

    class GravitationalModelParser {
    public:
        GravitationalModelParser(const std::string& filename);

        bool parse();

        int getMaxDegree() const;
        double getRadius() const;
        double getMu() const;

        // Returns the parsed data as a 2D vector of doubles
        const std::vector<std::vector<double>>& getData() const;

    private:
        std::string file;
        int Max_Degree;
        double Radius;
        double mu;
        std::vector<std::vector<double>> data;
    };
    
}