#include "common.h"
#include <string>

namespace orb {
InputParser::InputParser(const std::string& file) : file_(file) {
    parse();
}

std::string InputParser::trim(const std::string& str) {
    const std::string whitespace = " \t";
    size_t start = str.find_first_not_of(whitespace);
    if (start == std::string::npos) {
        return std::string(); // all whitespace
    }
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

} // namespace orb
