#include "common.h"

#include <string>
#include <unordered_map>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace orb {

ConfigParser::ConfigMap ConfigParser::parse(const std::string& text) {
    ConfigMap config;

    std::istringstream iss(text);
    std::string line;
    std::regex re(R"(^\s*([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\s*=\s*(.+?)\s*$)");

    while (std::getline(iss, line)) {
        if (line.empty()) continue; // skip empty lines

        std::smatch match;
        if (std::regex_match(line, match, re)) {
            std::string section = match[1];
            std::string key     = match[2];
            std::string valueStr = match[3];

            try {
                double value = std::stod(valueStr);
                config[section][key] = value;
            } catch (const std::invalid_argument&) {
                throw std::runtime_error(
                    "Invalid value for " + section + "." + key + ": " + valueStr
                );
            } catch (const std::out_of_range&) {
                throw std::runtime_error(
                    "Value out of range for " + section + "." + key + ": " + valueStr
                );
            }
        }
    }

    return config;
}

} // namespace orb
