#pragma once

// ExprProgram: the AST from Expr.h compiled to a flat postfix program.
//
// Recursive evalExpr() pays a call, a switch and a pointer chase per node per
// point. Grid sampling evaluates the same expression at (res+1)^3 points, so
// the extractor and the slice exporter instead compile the AST once and run
// it over a whole row of points at a time: every instruction becomes a tight
// loop over the row. Constant subtrees are folded at compile time, and
// constant integer exponents use exponentiation by squaring (as mc_pow does).
// The operations and their order are exactly those of evalExpr<double>, so
// the results are identical.

#include "Expr.h"

#include <algorithm>
#include <cstddef>
#include <vector>

class ExprProgram {
public:
    // One input coordinate of a batch: either a per-point array or a value
    // shared by every point of the batch.
    struct Lane {
        const double* values = nullptr;
        double scalar = 0.0;
        static Lane shared(double v) { return Lane{nullptr, v}; }
        static Lane array(const double* p) { return Lane{p, 0.0}; }
    };

    ExprProgram() = default;
    explicit ExprProgram(const Node* ast) {
        if (!ast) throw ExprError("ExprProgram: empty expression");
        int depth = 0;
        compile(ast, depth);
    }

    bool empty() const { return code_.empty(); }

    // Evaluates F at (x[i], y[i], z[i], k) for i < count. `scratch` is reused
    // between calls to avoid allocations.
    void evaluate(Lane x, Lane y, Lane z, double k, std::size_t count, double* out,
                  std::vector<double>& scratch) const {
        if (count == 0) return;
        scratch.resize(static_cast<std::size_t>(maxDepth_) * count);
        auto slot = [&](int s) { return scratch.data() + static_cast<std::size_t>(s) * count; };
        auto load = [&](double* dst, const Lane& lane) {
            if (lane.values) std::copy(lane.values, lane.values + count, dst);
            else std::fill(dst, dst + count, lane.scalar);
        };
        int sp = 0;
        for (const Instr& in : code_) {
            switch (in.op) {
                case Op::Const: std::fill(slot(sp), slot(sp) + count, in.value); ++sp; break;
                case Op::X: load(slot(sp++), x); break;
                case Op::Y: load(slot(sp++), y); break;
                case Op::Z: load(slot(sp++), z); break;
                case Op::K: std::fill(slot(sp), slot(sp) + count, k); ++sp; break;
                case Op::PowInt: {
                    double* a = slot(sp - 1);
                    for (std::size_t i = 0; i < count; ++i) a[i] = powInt(a[i], in.n);
                    break;
                }
                default:
                    if (isUnary(in.op)) {
                        double* a = slot(sp - 1);
                        for (std::size_t i = 0; i < count; ++i) a[i] = applyUnary(in.op, a[i]);
                    } else {
                        double* a = slot(sp - 2);
                        const double* b = slot(sp - 1);
                        binaryLoop(in.op, a, b, count);
                        --sp;
                    }
                    break;
            }
        }
        std::copy(slot(0), slot(0) + count, out);
    }

    // Single-point evaluation for any evalExpr-compatible number type.
    template <typename T>
    T evaluate(const T& x, const T& y, const T& z, const T& k) const {
        T stackBuffer[kInlineStack];
        std::vector<T> heap;
        T* stack = stackBuffer;
        if (maxDepth_ > kInlineStack) {
            heap.resize(static_cast<std::size_t>(maxDepth_));
            stack = heap.data();
        }
        int sp = 0;
        for (const Instr& in : code_) {
            switch (in.op) {
                case Op::Const: stack[sp++] = T(in.value); break;
                case Op::X: stack[sp++] = x; break;
                case Op::Y: stack[sp++] = y; break;
                case Op::Z: stack[sp++] = z; break;
                case Op::K: stack[sp++] = k; break;
                case Op::Neg: stack[sp - 1] = -stack[sp - 1]; break;
                case Op::Sin: stack[sp - 1] = mc_sin(stack[sp - 1]); break;
                case Op::Cos: stack[sp - 1] = mc_cos(stack[sp - 1]); break;
                case Op::Tan: stack[sp - 1] = mc_tan(stack[sp - 1]); break;
                case Op::Exp: stack[sp - 1] = mc_exp(stack[sp - 1]); break;
                case Op::Log: stack[sp - 1] = mc_log(stack[sp - 1]); break;
                case Op::Sqrt: stack[sp - 1] = mc_sqrt(stack[sp - 1]); break;
                case Op::Abs: stack[sp - 1] = mc_abs(stack[sp - 1]); break;
                case Op::PowInt: stack[sp - 1] = mc_pow(stack[sp - 1], T(static_cast<double>(in.n))); break;
                case Op::Add: --sp; stack[sp - 1] = stack[sp - 1] + stack[sp]; break;
                case Op::Sub: --sp; stack[sp - 1] = stack[sp - 1] - stack[sp]; break;
                case Op::Mul: --sp; stack[sp - 1] = stack[sp - 1] * stack[sp]; break;
                case Op::Div: --sp; stack[sp - 1] = stack[sp - 1] / stack[sp]; break;
                case Op::Pow: --sp; stack[sp - 1] = mc_pow(stack[sp - 1], stack[sp]); break;
                case Op::Min: --sp; stack[sp - 1] = mc_min(stack[sp - 1], stack[sp]); break;
                case Op::Max: --sp; stack[sp - 1] = mc_max(stack[sp - 1], stack[sp]); break;
            }
        }
        return stack[0];
    }

private:
    enum class Op : unsigned char {
        Const, X, Y, Z, K,
        Neg, Sin, Cos, Tan, Exp, Log, Sqrt, Abs, PowInt, // unary
        Add, Sub, Mul, Div, Pow, Min, Max                // binary
    };
    struct Instr {
        Op op;
        int n = 0;          // PowInt exponent
        double value = 0.0; // Const value
    };
    static constexpr int kInlineStack = 32;

    std::vector<Instr> code_;
    int maxDepth_ = 0;

    static bool isUnary(Op op) { return op >= Op::Neg && op <= Op::PowInt; }

    static bool hasVariable(const Node* n) {
        return exprContains(n, Node::VAR_X) || exprContains(n, Node::VAR_Y) ||
               exprContains(n, Node::VAR_Z) || exprContains(n, Node::VAR_K);
    }

    void emit(Instr in, int& depth, int delta) {
        code_.push_back(in);
        depth += delta;
        maxDepth_ = std::max(maxDepth_, depth);
    }

    void compile(const Node* n, int& depth) {
        if (n->type != Node::NUM && !hasVariable(n)) {
            // Constant subtree: same operations, evaluated once.
            emit({Op::Const, 0, evalExpr<double>(n, 0.0, 0.0, 0.0, 0.0)}, depth, +1);
            return;
        }
        switch (n->type) {
            case Node::NUM: emit({Op::Const, 0, n->num}, depth, +1); return;
            case Node::VAR_X: emit({Op::X}, depth, +1); return;
            case Node::VAR_Y: emit({Op::Y}, depth, +1); return;
            case Node::VAR_Z: emit({Op::Z}, depth, +1); return;
            case Node::VAR_K: emit({Op::K}, depth, +1); return;
            case Node::NEG: compile(n->kids[0].get(), depth); emit({Op::Neg}, depth, 0); return;
            case Node::ADD: compileBinary(n, Op::Add, depth); return;
            case Node::SUB: compileBinary(n, Op::Sub, depth); return;
            case Node::MUL: compileBinary(n, Op::Mul, depth); return;
            case Node::DIV: compileBinary(n, Op::Div, depth); return;
            case Node::POW: compilePow(n, depth); return;
            case Node::CALL:
                switch (n->func) {
                    case Node::SIN: compileUnary(n, Op::Sin, depth); return;
                    case Node::COS: compileUnary(n, Op::Cos, depth); return;
                    case Node::TAN: compileUnary(n, Op::Tan, depth); return;
                    case Node::EXP: compileUnary(n, Op::Exp, depth); return;
                    case Node::LOG: compileUnary(n, Op::Log, depth); return;
                    case Node::SQRT: compileUnary(n, Op::Sqrt, depth); return;
                    case Node::ABS: compileUnary(n, Op::Abs, depth); return;
                    case Node::FPOW: compilePow(n, depth); return;
                    case Node::MIN: compileBinary(n, Op::Min, depth); return;
                    case Node::MAX: compileBinary(n, Op::Max, depth); return;
                    case Node::NONE: break;
                }
                break;
        }
        throw ExprError("internal: unhandled node type");
    }

    void compileUnary(const Node* n, Op op, int& depth) {
        compile(n->kids[0].get(), depth);
        emit({op}, depth, 0);
    }
    void compileBinary(const Node* n, Op op, int& depth) {
        compile(n->kids[0].get(), depth);
        compile(n->kids[1].get(), depth);
        emit({op}, depth, -1);
    }
    void compilePow(const Node* n, int& depth) {
        const Node* exponent = n->kids[1].get();
        int power = 0;
        if (!hasVariable(exponent) &&
            smallIntegerExponent(evalExpr<double>(exponent, 0.0, 0.0, 0.0, 0.0), power)) {
            compile(n->kids[0].get(), depth);
            emit({Op::PowInt, power}, depth, 0);
            return;
        }
        compileBinary(n, Op::Pow, depth);
    }

    static double applyUnary(Op op, double a) {
        switch (op) {
            case Op::Neg: return -a;
            case Op::Sin: return mc_sin(a);
            case Op::Cos: return mc_cos(a);
            case Op::Tan: return mc_tan(a);
            case Op::Exp: return mc_exp(a);
            case Op::Log: return mc_log(a);
            case Op::Sqrt: return mc_sqrt(a);
            case Op::Abs: return mc_abs(a);
            default: return a;
        }
    }

    static void binaryLoop(Op op, double* a, const double* b, std::size_t count) {
        switch (op) {
            case Op::Add: for (std::size_t i = 0; i < count; ++i) a[i] = a[i] + b[i]; break;
            case Op::Sub: for (std::size_t i = 0; i < count; ++i) a[i] = a[i] - b[i]; break;
            case Op::Mul: for (std::size_t i = 0; i < count; ++i) a[i] = a[i] * b[i]; break;
            case Op::Div: for (std::size_t i = 0; i < count; ++i) a[i] = a[i] / b[i]; break;
            case Op::Pow: for (std::size_t i = 0; i < count; ++i) a[i] = mc_pow(a[i], b[i]); break;
            case Op::Min: for (std::size_t i = 0; i < count; ++i) a[i] = mc_min(a[i], b[i]); break;
            case Op::Max: for (std::size_t i = 0; i < count; ++i) a[i] = mc_max(a[i], b[i]); break;
            default: break;
        }
    }
};
