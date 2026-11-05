#pragma once

// File-system-safe names derived from user input.

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>

namespace iso_files {

// Converts an expression into a portable, readable file or folder name
// component ('*' becomes 'x', whitespace and unsafe characters '_').
// Callers keep the original expression in the file contents.
inline std::string safeFileComponent(const std::string& expression, std::size_t maxLength = 64) {
    std::string result;
    result.reserve(std::min(expression.size(), maxLength));
    bool previousUnderscore = false;
    for (unsigned char c : expression) {
        char out = '_';
        if (std::isalnum(c) || c == '_' || c == '-' || c == '+' ||
            c == '^' || c == '(' || c == ')' || c == '.') {
            out = static_cast<char>(c);
        } else if (c == '*') {
            out = 'x'; // multiplication remains visually distinct in folder names
        } else if (std::isspace(c)) {
            out = '_';
        }
        if (out == '_') {
            if (previousUnderscore) continue;
            previousUnderscore = true;
        } else {
            previousUnderscore = false;
        }
        result.push_back(out);
        if (result.size() >= maxLength) break;
    }
    while (!result.empty() && (result.back() == '.' || result.back() == ' '))
        result.pop_back();
    while (!result.empty() && result.front() == '.')
        result.erase(result.begin());
    if (result.empty()) result = "expression";
    if (result == "CON" || result == "PRN" || result == "AUX" ||
        result == "NUL" || result == "COM1" || result == "LPT1") {
        result = "expr_" + result;
    }
    return result;
}

} // namespace iso_files
