#pragma once

// Text formatting shared by every CSV/report writer.

#include <charconv>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>

namespace iso_csv {

// Fast appenders for large exports. append(out, double) produces exactly what
// an ostream with setprecision(17) prints ("%.17g"), without stream overhead.
// std::to_chars(general, 17) is specified to match printf("%.17g") and is
// several times faster; CMake defines ISO_HAVE_FLOAT_TO_CHARS when the
// standard library provides it (libc++ does not advertise it via
// __cpp_lib_to_chars), otherwise snprintf is used.
inline void append(std::string& out, double value) {
    char buffer[40];
#if defined(ISO_HAVE_FLOAT_TO_CHARS) || (defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L)
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::general, 17);
    out.append(buffer, result.ptr);
#else
    const int n = std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    if (n > 0) out.append(buffer, static_cast<std::size_t>(n));
#endif
}
template <typename Int, typename = std::enable_if_t<std::is_integral_v<Int>>>
inline void append(std::string& out, Int value) {
    char buffer[24];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    out.append(buffer, result.ptr);
}
inline void append(std::string& out, const char* text) { out += text; }

// RFC 4180 quoting: wrap in double quotes and double any embedded quote.
inline std::string quote(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 2);
    out += '"';
    for (char ch : value) {
        if (ch == '"') out += '"';
        out += ch;
    }
    out += '"';
    return out;
}

inline std::string formatWithPrecision(double value, int digits) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream out;
    out << std::setprecision(digits) << value;
    return out.str();
}

// Round-trip precision for IEEE-754 doubles; non-finite values become "NA".
inline std::string number(double value) { return formatWithPrecision(value, 17); }

// Compact form for on-screen tables.
inline std::string shortNumber(double value) { return formatWithPrecision(value, 8); }

} // namespace iso_csv
