#include "analysis/SingularityBatch.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}
}

int main() {
    const std::string gLines =
        "1\nx\ny\nz\nx+y\nx+z\ny+z\nx-y\nx-z\ny-z\n"
        "x^2\ny^2\nz^2\nx*y\nx*z\ny*z\n"
        "x^2+y^2+z^2\nx*y+y*z+z*x\nx^3\ny^3\nz^3\n"
        "x*y*z\nx^2*y\nx^2*z\ny^2*z\n";

    const auto expressions = iso_batch::parseGExpressions(gLines);
    require(expressions.size() == 25, "exactly 25 G expressions are parsed");
    require(expressions.front() == "1", "first G expression is preserved");
    require(expressions.back() == "y^2*z", "last G expression is preserved");

    bool rejectedCount = false;
    try {
        (void)iso_batch::parseGExpressions("1\nx\n");
    } catch (const std::invalid_argument&) {
        rejectedCount = true;
    }
    require(rejectedCount, "a list other than 25 G expressions is rejected");

    const std::string family = iso_batch::composeFamilyExpression(
        "x^3+y^3+z^3", "x*y");
    NodePtr ast = parseExpression(family);
    const auto affine = iso_convergence::decomposeAffineK(ast.get());
    require(affine.valid, "constructed F+kG expression is affine in k");
    require(iso_batch::containsKVariable(ast.get()),
            "combined family contains k as the intended parameter");

    NodePtr base = parseExpression("x^3+y^3+z^3");
    NodePtr perturbation = parseExpression("x*y");
    require(!iso_batch::containsKVariable(base.get()), "F without k is accepted");
    require(!iso_batch::containsKVariable(perturbation.get()), "G without k is accepted");
    NodePtr invalid = parseExpression("x+k");
    require(iso_batch::containsKVariable(invalid.get()), "k inside an input component is detected");

    const std::string safe = iso_batch::safeFolderComponent("x/y + z*z");
    require(safe.find('/') == std::string::npos, "folder names cannot contain path separators");
    require(!safe.empty(), "folder names are non-empty");

    std::cout << "Batch helper tests passed.\n";
    return 0;
}
