// Expr.h
//
// A minimal recursive-descent parser + AST for scalar expressions in the
// variables x, y, z, k. The AST is evaluated with a template so the exact
// same tree works for T = double (fast path, just a value), T = Dual3
// (value + exact gradient via forward-mode AD, see Dual3.h) and the
// third-order Taylor jets used by SingularityResearch.h.
//
// Function names are resolved to Node::Func once, at parse time, so
// evaluation never compares strings and an unknown function or a wrong
// argument count is reported as a parse error instead of failing later in
// the middle of a mesh build.
#pragma once
#include "Dual3.h"
#include <cctype>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct ExprError : std::runtime_error {
    explicit ExprError(const std::string& msg) : std::runtime_error(msg) {}
};

struct Node {
    enum Type { NUM, VAR_X, VAR_Y, VAR_Z, VAR_K, NEG, ADD, SUB, MUL, DIV, POW, CALL };
    enum Func : std::uint8_t { NONE, SIN, COS, TAN, EXP, LOG, SQRT, ABS, FPOW, MIN, MAX };
    Type type = NUM;
    Func func = NONE;      // only for CALL
    double num = 0.0;      // only for NUM
    std::string fname;     // only for CALL: lower-cased name as typed (e.g. "ln"), used for text output
    std::vector<std::unique_ptr<Node>> kids;
};

using NodePtr = std::unique_ptr<Node>;

namespace exprimpl {

constexpr double kPi = 3.14159265358979323846;
constexpr double kE = 2.71828182845904523536;

struct FunctionInfo {
    const char* name;
    Node::Func func;
    std::size_t arity;
};

inline const FunctionInfo* lookupFunction(const std::string& name) {
    static const FunctionInfo table[] = {
        {"sin", Node::SIN, 1},  {"cos", Node::COS, 1},   {"tan", Node::TAN, 1},
        {"exp", Node::EXP, 1},  {"log", Node::LOG, 1},   {"ln", Node::LOG, 1},
        {"sqrt", Node::SQRT, 1}, {"abs", Node::ABS, 1},  {"pow", Node::FPOW, 2},
        {"min", Node::MIN, 2},  {"max", Node::MAX, 2},
    };
    for (const FunctionInfo& info : table)
        if (name == info.name) return &info;
    return nullptr;
}

struct Token {
    enum Kind { END, NUMBER, IDENT, PLUS, MINUS, STAR, SLASH, CARET, LPAREN, RPAREN, COMMA };
    Kind kind = END;
    double num = 0.0;
    std::string text;
};

class Lexer {
public:
    explicit Lexer(const std::string& s) : src(s) { advance(); }
    const Token& peek() const { return cur; }
    Token take() { Token t = std::move(cur); advance(); return t; }

private:
    const std::string& src;
    std::size_t pos = 0;
    Token cur;

    bool digitAt(std::size_t i) const {
        return i < src.size() && std::isdigit(static_cast<unsigned char>(src[i]));
    }

    void advance() {
        while (pos < src.size() && std::isspace(static_cast<unsigned char>(src[pos]))) ++pos;
        cur = Token{};
        if (pos >= src.size()) return;
        const char c = src[pos];
        if (digitAt(pos) || (c == '.' && digitAt(pos + 1))) {
            const std::size_t start = pos;
            bool sawDot = false, sawExp = false;
            while (pos < src.size()) {
                const char d = src[pos];
                if (std::isdigit(static_cast<unsigned char>(d))) { ++pos; }
                else if (d == '.' && !sawDot && !sawExp) { sawDot = true; ++pos; }
                else if ((d == 'e' || d == 'E') && !sawExp &&
                         (digitAt(pos + 1) ||
                          ((pos + 1 < src.size() && (src[pos + 1] == '+' || src[pos + 1] == '-')) && digitAt(pos + 2)))) {
                    sawExp = true;
                    ++pos;
                    if (src[pos] == '+' || src[pos] == '-') ++pos;
                } else break;
            }
            cur.kind = Token::NUMBER;
            cur.num = std::stod(src.substr(start, pos - start));
            return;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            const std::size_t start = pos;
            while (pos < src.size() && (std::isalnum(static_cast<unsigned char>(src[pos])) || src[pos] == '_')) ++pos;
            cur.kind = Token::IDENT;
            cur.text = src.substr(start, pos - start);
            return;
        }
        ++pos;
        switch (c) {
            case '+': cur.kind = Token::PLUS; return;
            case '-': cur.kind = Token::MINUS; return;
            case '*': cur.kind = Token::STAR; return;
            case '/': cur.kind = Token::SLASH; return;
            case '^': cur.kind = Token::CARET; return;
            case '(': cur.kind = Token::LPAREN; return;
            case ')': cur.kind = Token::RPAREN; return;
            case ',': cur.kind = Token::COMMA; return;
            default:
                throw ExprError(std::string("unexpected character '") + c + "'");
        }
    }
};

class Parser {
public:
    explicit Parser(const std::string& s) : lex(s) {}

    NodePtr parseTop() {
        NodePtr n = parseExpr();
        if (lex.peek().kind != Token::END)
            throw ExprError("unexpected trailing input");
        return n;
    }

private:
    Lexer lex;

    static NodePtr makeNode(Node::Type t) {
        auto n = std::make_unique<Node>();
        n->type = t;
        return n;
    }

    static NodePtr makeBin(Node::Type t, NodePtr a, NodePtr b) {
        NodePtr n = makeNode(t);
        n->kids.push_back(std::move(a));
        n->kids.push_back(std::move(b));
        return n;
    }

    bool accept(Token::Kind kind) {
        if (lex.peek().kind != kind) return false;
        lex.take();
        return true;
    }

    // expr := term (('+'|'-') term)*
    NodePtr parseExpr() {
        NodePtr lhs = parseTerm();
        for (;;) {
            if (accept(Token::PLUS)) lhs = makeBin(Node::ADD, std::move(lhs), parseTerm());
            else if (accept(Token::MINUS)) lhs = makeBin(Node::SUB, std::move(lhs), parseTerm());
            else return lhs;
        }
    }

    // term := power (('*'|'/') power)*
    NodePtr parseTerm() {
        NodePtr lhs = parsePower();
        for (;;) {
            if (accept(Token::STAR)) lhs = makeBin(Node::MUL, std::move(lhs), parsePower());
            else if (accept(Token::SLASH)) lhs = makeBin(Node::DIV, std::move(lhs), parsePower());
            else return lhs;
        }
    }

    // power := unary ('^' power)?     (right-associative)
    NodePtr parsePower() {
        NodePtr base = parseUnary();
        if (accept(Token::CARET)) return makeBin(Node::POW, std::move(base), parsePower());
        return base;
    }

    // unary := '-' unary | '+' unary | primary
    NodePtr parseUnary() {
        if (accept(Token::MINUS)) {
            NodePtr n = makeNode(Node::NEG);
            n->kids.push_back(parseUnary());
            return n;
        }
        if (accept(Token::PLUS)) return parseUnary();
        return parsePrimary();
    }

    // primary := NUMBER | IDENT | IDENT '(' args ')' | '(' expr ')'
    NodePtr parsePrimary() {
        const Token::Kind kind = lex.peek().kind;
        if (kind == Token::NUMBER) {
            NodePtr n = makeNode(Node::NUM);
            n->num = lex.take().num;
            return n;
        }
        if (kind == Token::LPAREN) {
            lex.take();
            NodePtr n = parseExpr();
            if (!accept(Token::RPAREN)) throw ExprError("expected ')'");
            return n;
        }
        if (kind == Token::IDENT) {
            const std::string typed = lex.take().text;
            std::string name = typed;
            for (char& c : name) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (accept(Token::LPAREN)) return parseCall(name);

            NodePtr n = makeNode(Node::NUM);
            if (name == "x") n->type = Node::VAR_X;
            else if (name == "y") n->type = Node::VAR_Y;
            else if (name == "z") n->type = Node::VAR_Z;
            else if (name == "k") n->type = Node::VAR_K;
            else if (name == "pi") n->num = kPi;
            else if (name == "e") n->num = kE;
            else throw ExprError("unknown identifier '" + typed + "' (use x, y, z, k)");
            return n;
        }
        throw ExprError("expected a number, variable or '('");
    }

    NodePtr parseCall(const std::string& name) {
        NodePtr n = makeNode(Node::CALL);
        n->fname = name;
        if (lex.peek().kind != Token::RPAREN) {
            n->kids.push_back(parseExpr());
            while (accept(Token::COMMA)) n->kids.push_back(parseExpr());
        }
        if (!accept(Token::RPAREN)) throw ExprError("expected ')' after function arguments");

        const FunctionInfo* info = lookupFunction(name);
        if (!info)
            throw ExprError("unknown function '" + name + "'");
        if (n->kids.size() != info->arity)
            throw ExprError("function '" + name + "' expects " + std::to_string(info->arity) +
                            (info->arity == 1 ? " argument" : " arguments") + ", got " +
                            std::to_string(n->kids.size()));
        n->func = info->func;
        return n;
    }
};

} // namespace exprimpl

// Parse `src` into an AST. Throws ExprError on malformed input.
inline NodePtr parseExpression(const std::string& src) {
    exprimpl::Parser p(src);
    return p.parseTop();
}

// True when the subtree references the given variable node type.
inline bool exprContains(const Node* n, Node::Type variable) {
    if (!n) return false;
    if (n->type == variable) return true;
    for (const auto& child : n->kids)
        if (exprContains(child.get(), variable)) return true;
    return false;
}

// Evaluate the AST at (x,y,z,k). T is either `double` (value only) or
// `Dual3` (value + exact partial derivatives, when x/y/z are seeded as
// Dual3(value, 1,0,0 / 0,1,0 / 0,0,1) and k as Dual3(value) with zero
// derivative). Any other T works if it provides the arithmetic operators
// and the mc_* elementary functions.
template <typename T>
T evalExpr(const Node* n, const T& x, const T& y, const T& z, const T& k) {
    auto kid = [&](std::size_t i) { return evalExpr<T>(n->kids[i].get(), x, y, z, k); };
    switch (n->type) {
        case Node::NUM:   return T(n->num);
        case Node::VAR_X: return x;
        case Node::VAR_Y: return y;
        case Node::VAR_Z: return z;
        case Node::VAR_K: return k;
        case Node::NEG:   return -kid(0);
        case Node::ADD:   return kid(0) + kid(1);
        case Node::SUB:   return kid(0) - kid(1);
        case Node::MUL:   return kid(0) * kid(1);
        case Node::DIV:   return kid(0) / kid(1);
        case Node::POW:   return mc_pow(kid(0), kid(1));
        case Node::CALL:
            switch (n->func) {
                case Node::SIN:  return mc_sin(kid(0));
                case Node::COS:  return mc_cos(kid(0));
                case Node::TAN:  return mc_tan(kid(0));
                case Node::EXP:  return mc_exp(kid(0));
                case Node::LOG:  return mc_log(kid(0));
                case Node::SQRT: return mc_sqrt(kid(0));
                case Node::ABS:  return mc_abs(kid(0));
                case Node::FPOW: return mc_pow(kid(0), kid(1));
                case Node::MIN:  return mc_min(kid(0), kid(1));
                case Node::MAX:  return mc_max(kid(0), kid(1));
                case Node::NONE: break;
            }
            break;
    }
    throw ExprError("internal: unhandled node type");
}
