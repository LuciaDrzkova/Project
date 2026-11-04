#pragma once

// Runs a batch of families F(x,y,z) + k G_i(x,y,z), exporting the same four
// CSV reports as the interactive single-run convergence action. Families run
// one after another; within a family every k sample is processed in parallel.

#include "SingularityResearch.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iso_batch {

inline constexpr std::size_t kRequiredGCount = 25;

struct Settings {
    std::string fExpression;
    std::vector<std::string> gExpressions;
    double kMinimum = -0.1;
    double kMaximum = 0.1;
    int kSampleCount = 25;
    int resolution = 48;
    double halfExtent = 1.0;
    iso_convergence::Point3 probe{};
};

struct RunReport {
    std::string gExpression;
    std::filesystem::path directory;
    bool success = false;
    std::string message;
};

struct BatchResult {
    std::filesystem::path fDirectory;
    std::vector<RunReport> runs;

    std::size_t successCount() const {
        return static_cast<std::size_t>(std::count_if(
            runs.begin(), runs.end(), [](const RunReport& run) { return run.success; }));
    }

    std::size_t failureCount() const { return runs.size() - successCount(); }
};

inline std::string trim(const std::string& value) {
    const auto first = std::find_if_not(value.begin(), value.end(),
        [](unsigned char c) { return std::isspace(c) != 0; });
    if (first == value.end()) return {};
    const auto last = std::find_if_not(value.rbegin(), value.rend(),
        [](unsigned char c) { return std::isspace(c) != 0; }).base();
    return std::string(first, last);
}

// One non-empty, non-comment line corresponds to one G expression.
inline std::vector<std::string> parseGExpressions(const std::string& multiline) {
    std::vector<std::string> expressions;
    std::istringstream input(multiline);
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty() || line.front() == '#') continue;
        expressions.push_back(std::move(line));
    }
    if (expressions.size() != kRequiredGCount) {
        throw std::invalid_argument(
            "Enter exactly 25 non-empty G expressions, one expression per line. "
            "Lines beginning with # are ignored. Found " +
            std::to_string(expressions.size()) + ".");
    }
    return expressions;
}

inline std::string composeFamilyExpression(const std::string& fExpression,
                                           const std::string& gExpression) {
    return "(" + trim(fExpression) + ")+k*(" + trim(gExpression) + ")";
}

inline bool containsKVariable(const Node* node) { return exprContains(node, Node::VAR_K); }

// Convert an expression into a portable, readable folder component. The
// original unsanitized expressions remain recorded in the CSV metadata.
inline std::string safeFolderComponent(const std::string& expression,
                                       std::size_t maxLength = 64) {
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

// Readable, expression-based run folder; name collisions get _02, _03, ... suffixes.
inline std::filesystem::path uniqueRunDirectory(
    const std::filesystem::path& parent,
    const std::string& fExpression,
    const std::string& gExpression) {
    const std::string fPart = safeFolderComponent(fExpression, 34);
    const std::string gPart = safeFolderComponent(gExpression, 34);
    const std::string base = "run_" + fPart + "_plus_k_" + gPart;
    std::filesystem::path candidate = parent / base;
    std::size_t suffix = 2;
    while (std::filesystem::exists(candidate)) {
        std::ostringstream name;
        name << base << "_" << std::setw(2) << std::setfill('0') << suffix++;
        candidate = parent / name.str();
    }
    return candidate;
}

inline void validateSettings(const Settings& settings) {
    if (trim(settings.fExpression).empty())
        throw std::invalid_argument("F(x,y,z) must not be empty.");
    if (settings.gExpressions.size() != kRequiredGCount)
        throw std::invalid_argument("Exactly 25 G expressions are required.");
    if (!std::isfinite(settings.kMinimum) || !std::isfinite(settings.kMaximum) ||
        !(settings.kMinimum < settings.kMaximum))
        throw std::invalid_argument("k minimum and maximum must be finite, with minimum < maximum.");
    if (settings.kSampleCount < 3 || settings.kSampleCount > 101)
        throw std::invalid_argument("Number of k samples must be between 3 and 101.");
    if (settings.resolution < 16 || settings.resolution > 110)
        throw std::invalid_argument("Mesh resolution must be between 16 and 110 cells per axis.");
    if (!(settings.halfExtent > 0.0) || !std::isfinite(settings.halfExtent))
        throw std::invalid_argument("Domain half-extent must be finite and positive.");
    if (!std::isfinite(settings.probe.x) || !std::isfinite(settings.probe.y) ||
        !std::isfinite(settings.probe.z))
        throw std::invalid_argument("Point of interest P must contain finite coordinates.");
    for (std::size_t i = 0; i < settings.gExpressions.size(); ++i) {
        if (trim(settings.gExpressions[i]).empty())
            throw std::invalid_argument("G expression " + std::to_string(i + 1) + " is empty.");
    }
}

inline BatchResult runBatch(const Settings& settings,
                            const std::filesystem::path& experimentsDirectory,
                            std::atomic<int>* completedCounter = nullptr) {
    validateSettings(settings);

    // Validate F once, then each G separately. This prevents a hidden k in F or
    // G from changing the intended meaning of F + kG.
    NodePtr fAst = parseExpression(settings.fExpression);
    if (containsKVariable(fAst.get()))
        throw std::invalid_argument("F must depend only on x, y and z; do not include k in F.");

    BatchResult result;
    result.fDirectory = experimentsDirectory / safeFolderComponent(settings.fExpression, 80);
    std::filesystem::create_directories(result.fDirectory);
    result.runs.reserve(settings.gExpressions.size());

    const std::vector<double> kValues = iso_convergence::makeLinearKValues(
        settings.kMinimum, settings.kMaximum, settings.kSampleCount);

    for (std::size_t i = 0; i < settings.gExpressions.size(); ++i) {
        RunReport report;
        report.gExpression = settings.gExpressions[i];
        report.directory = uniqueRunDirectory(result.fDirectory, settings.fExpression, report.gExpression);

        try {
            NodePtr gAst = parseExpression(report.gExpression);
            if (containsKVariable(gAst.get()))
                throw std::invalid_argument("G must depend only on x, y and z; do not include k in G.");

            const std::string familyExpression = composeFamilyExpression(
                settings.fExpression, report.gExpression);
            NodePtr familyAst = parseExpression(familyExpression);
            if (!iso_convergence::decomposeAffineK(familyAst.get()).valid)
                throw std::invalid_argument("The constructed family could not be verified as affine in k.");

            iso_convergence::Options options;
            options.probe = settings.probe;
            options.resolution = settings.resolution;
            options.halfExtent = settings.halfExtent;
            options.localRadius = 0.45 * settings.halfExtent;
            options.kValues = kValues;
            options.familyExpression = familyExpression;

            // One shared, parallel mesh pass feeds both analyses.
            const iso_research::FamilyAnalysis analysis =
                iso_research::analyzeFamily(familyAst.get(), options, 3);
            const iso_convergence::Analysis& convergence = analysis.convergence;
            const iso_research::ResearchAnalysis& research = analysis.research;

            std::filesystem::create_directories(report.directory);
            std::string writeError;
            const auto filePath = [&report](const char* name) {
                return (report.directory / name).string();
            };
            if (!iso_convergence::writeCsv(
                    filePath("singularity_convergence.csv"), convergence, &writeError))
                throw std::runtime_error(writeError.empty() ? "Could not write convergence CSV." : writeError);
            if (!iso_research::writeTaylorCsv(
                    filePath("singularity_taylor.csv"), research, &writeError))
                throw std::runtime_error(writeError.empty() ? "Could not write Taylor CSV." : writeError);
            if (!iso_research::writeCriticalPointsCsv(
                    filePath("singularity_critical_points.csv"), research, &writeError))
                throw std::runtime_error(writeError.empty() ? "Could not write critical-points CSV." : writeError);
            if (!iso_research::writeTriangleGeometryCsv(
                    filePath("singularity_triangle_geometry.csv"), research, &writeError))
                throw std::runtime_error(writeError.empty() ? "Could not write triangle-geometry CSV." : writeError);

            report.success = true;
            report.message = research.summary;
        } catch (const std::exception& error) {
            report.success = false;
            report.message = error.what();
            std::error_code cleanupError;
            std::filesystem::remove_all(report.directory, cleanupError);
        } catch (...) {
            report.success = false;
            report.message = "Unknown error while analysing this F+kG family.";
            std::error_code cleanupError;
            std::filesystem::remove_all(report.directory, cleanupError);
        }

        result.runs.push_back(std::move(report));
        if (completedCounter)
            completedCounter->fetch_add(1, std::memory_order_relaxed);
    }

    return result;
}

inline std::string summarize(const BatchResult& result) {
    std::ostringstream out;
    out << "Batch finished: " << result.successCount() << " of " << result.runs.size()
        << " runs saved successfully.\nOutput: " << result.fDirectory.string();
    for (std::size_t i = 0; i < result.runs.size(); ++i) {
        if (!result.runs[i].success) {
            out << "\nG" << (i + 1) << " failed (" << result.runs[i].gExpression
                << "): " << result.runs[i].message;
        }
    }
    return out.str();
}

} // namespace iso_batch
