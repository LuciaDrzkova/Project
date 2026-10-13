// Expr.h
//
// A minimal recursive-descent parser + AST for scalar expressions in the
// variables x, y, z, k. The AST is evaluated with a template so the exact
// same tree works for T = double (fast path, just a value) and
// T = Dual3 (value + exact gradient via forward-mode AD, see Dual3.h).
#pragma once
#include "Dual3.h"
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cctype>
#include <cmath>

struct ExprError : std::runtime_error {
    explicit ExprError(const std::string& msg) : std::runtime_error(msg) {}
};

struct Node {
    enum Type { NUM, VAR_X, VAR_Y, VAR_Z, VAR_K, NEG, ADD, SUB, MUL, DIV, POW, CALL };
    Type type;
    double num = 0.0;
    std::string fname; // only for CALL
    std::vector<std::unique_ptr<Node>> kids;
};

using NodePtr = std::unique_ptr<Node>;

namespace exprimpl {

struct Token {
    enum Kind { END, NUMBER, IDENT, PLUS, MINUS, STAR, SLASH, CARET, LPAREN, RPAREN, COMMA };
    Kind kind;
    double num = 0.0;
    std::string text;
};

class Lexer {
public:
    explicit Lexer(const std::string& s) : src(s), pos(0) { advance(); }
    const Token& peek() const { return cur; }
    Token take() { Token t = cur; advance(); return t; }

private:
    std::string src;
    size_t pos;
    Token cur;

    void advance() {
        while (pos < src.size() && std::isspace((unsigned char)src[pos])) pos++;
        if (pos >= src.size()) { cur = { Token::END }; return; }
        char c = src[pos];
        if (std::isdigit((unsigned char)c) || (c == '.' && pos + 1 < src.size() && std::isdigit((unsigned char)src[pos+1]))) {
            size_t start = pos;
            bool sawDot = false, sawExp = false;
            while (pos < src.size()) {
                char d = src[pos];
                if (std::isdigit((unsigned char)d)) { pos++; }
                else if (d == '.' && !sawDot && !sawExp) { sawDot = true; pos++; }
                else if ((d == 'e' || d == 'E') && !sawExp &&
                         pos + 1 < src.size() &&
                         (std::isdigit((unsigned char)src[pos+1]) ||
                          ((src[pos+1]=='+'||src[pos+1]=='-') && pos+2<src.size() && std::isdigit((unsigned char)src[pos+2])))) {
                    sawExp = true; pos++;
                    if (src[pos]=='+'||src[pos]=='-') pos++;
                } else break;
            }
            std::string numtext = src.substr(start, pos - start);
            cur = { Token::NUMBER };
            cur.num = std::stod(numtext);
            return;
        }
        if (std::isalpha((unsigned char)c) || c == '_') {
            size_t start = pos;
            while (pos < src.size() && (std::isalnum((unsigned char)src[pos]) || src[pos] == '_')) pos++;
            cur = { Token::IDENT };
            cur.text = src.substr(start, pos - start);
            return;
        }
        pos++;
        switch (c) {
            case '+': cur = { Token::PLUS }; return;
            case '-': cur = { Token::MINUS }; return;
            case '*': cur = { Token::STAR }; return;
            case '/': cur = { Token::SLASH }; return;
            case '^': cur = { Token::CARET }; return;
            case '(': cur = { Token::LPAREN }; return;
            case ')': cur = { Token::RPAREN }; return;
            case ',': cur = { Token::COMMA }; return;
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

    static NodePtr makeBin(Node::Type t, NodePtr a, NodePtr b) {
        auto n = std::make_unique<Node>();
        n->type = t;
        n->kids.push_back(std::move(a));
        n->kids.push_back(std::move(b));
        return n;
    }

    // expr := term (('+'|'-') term)*
    NodePtr parseExpr() {
        NodePtr lhs = parseTerm();
        for (;;) {
            if (lex.peek().kind == Token::PLUS) { lex.take(); lhs = makeBin(Node::ADD, std::move(lhs), parseTerm()); }
            else if (lex.peek().kind == Token::MINUS) { lex.take(); lhs = makeBin(Node::SUB, std::move(lhs), parseTerm()); }
            else break;
        }
        return lhs;
    }

    // term := power (('*'|'/') power)*
    NodePtr parseTerm() {
        NodePtr lhs = parsePower();
        for (;;) {
            if (lex.peek().kind == Token::STAR) { lex.take(); lhs = makeBin(Node::MUL, std::move(lhs), parsePower()); }
            else if (lex.peek().kind == Token::SLASH) { lex.take(); lhs = makeBin(Node::DIV, std::move(lhs), parsePower()); }
            else break;
        }
        return lhs;
    }

    // power := unary ('^' power)?     (right-associative)
    NodePtr parsePower() {
        NodePtr base = parseUnary();
        if (lex.peek().kind == Token::CARET) {
            lex.take();
            NodePtr exp = parsePower();
            return makeBin(Node::POW, std::move(base), std::move(exp));
        }
        return base;
    }

    // unary := '-' unary | '+' unary | primary
    NodePtr parseUnary() {
        if (lex.peek().kind == Token::MINUS) {
            lex.take();
            auto n = std::make_unique<Node>();
            n->type = Node::NEG;
            n->kids.push_back(parseUnary());
            return n;
        }
        if (lex.peek().kind == Token::PLUS) {
            lex.take();
            return parseUnary();
        }
        return parsePrimary();
    }

    // primary := NUMBER | IDENT | IDENT '(' args ')' | '(' expr ')'
    NodePtr parsePrimary() {
        Token t = lex.peek();
        if (t.kind == Token::NUMBER) {
            lex.take();
            auto n = std::make_unique<Node>();
            n->type = Node::NUM;
            n->num = t.num;
            return n;
        }
        if (t.kind == Token::LPAREN) {
            lex.take();
            NodePtr n = parseExpr();
            if (lex.peek().kind != Token::RPAREN) throw ExprError("expected ')'");
            lex.take();
            return n;
        }
        if (t.kind == Token::IDENT) {
            lex.take();
            std::string name = t.text;
            for (auto& c : name) c = (char)std::tolower((unsigned char)c);
            if (lex.peek().kind == Token::LPAREN) {
                lex.take();
                auto n = std::make_unique<Node>();
                n->type = Node::CALL;
                n->fname = name;
                if (lex.peek().kind != Token::RPAREN) {
                    n->kids.push_back(parseExpr());
                    while (lex.peek().kind == Token::COMMA) {
                        lex.take();
                        n->kids.push_back(parseExpr());
                    }
                }
                if (lex.peek().kind != Token::RPAREN) throw ExprError("expected ')' after function arguments");
                lex.take();
                return n;
            }
            auto n = std::make_unique<Node>();
            if (name == "x") n->type = Node::VAR_X;
            else if (name == "y") n->type = Node::VAR_Y;
            else if (name == "z") n->type = Node::VAR_Z;
            else if (name == "k") n->type = Node::VAR_K;
            else if (name == "pi") { n->type = Node::NUM; n->num = M_PI; }
            else if (name == "e")  { n->type = Node::NUM; n->num = M_E; }
            else throw ExprError("unknown identifier '" + t.text + "' (use x, y, z, k)");
            return n;
        }
        throw ExprError("expected a number, variable or '('");
    }
};

} // namespace exprimpl

// Parse `src` into an AST. Throws ExprError on malformed input.
inline NodePtr parseExpression(const std::string& src) {
    exprimpl::Parser p(src);
    return p.parseTop();
}

// Evaluate the AST at (x,y,z,k). T is either `double` (value only) or
// `Dual3` (value + exact partial derivatives, when x/y/z are seeded as
// Dual3(value, 1,0,0 / 0,1,0 / 0,0,1) and k as Dual3(value) with zero
// derivative).
template <typename T>
T evalExpr(const Node* n, const T& x, const T& y, const T& z, const T& k) {
    switch (n->type) {
        case Node::NUM:   return T(n->num);
        case Node::VAR_X: return x;
        case Node::VAR_Y: return y;
        case Node::VAR_Z: return z;
        case Node::VAR_K: return k;
        case Node::NEG:   return -evalExpr<T>(n->kids[0].get(), x, y, z, k);
        case Node::ADD:   return evalExpr<T>(n->kids[0].get(), x, y, z, k) + evalExpr<T>(n->kids[1].get(), x, y, z, k);
        case Node::SUB:   return evalExpr<T>(n->kids[0].get(), x, y, z, k) - evalExpr<T>(n->kids[1].get(), x, y, z, k);
        case Node::MUL:   return evalExpr<T>(n->kids[0].get(), x, y, z, k) * evalExpr<T>(n->kids[1].get(), x, y, z, k);
        case Node::DIV:   return evalExpr<T>(n->kids[0].get(), x, y, z, k) / evalExpr<T>(n->kids[1].get(), x, y, z, k);
        case Node::POW:   return mc_pow(evalExpr<T>(n->kids[0].get(), x, y, z, k), evalExpr<T>(n->kids[1].get(), x, y, z, k));
        case Node::CALL: {
            const auto& f = n->fname;
            if (n->kids.size() == 1) {
                T a = evalExpr<T>(n->kids[0].get(), x, y, z, k);
                if (f == "sin")  return mc_sin(a);
                if (f == "cos")  return mc_cos(a);
                if (f == "tan")  return mc_tan(a);
                if (f == "exp")  return mc_exp(a);
                if (f == "log" || f == "ln") return mc_log(a);
                if (f == "sqrt") return mc_sqrt(a);
                if (f == "abs")  return mc_abs(a);
                throw ExprError("unknown function '" + f + "' (1 argument)");
            }
            if (n->kids.size() == 2) {
                T a = evalExpr<T>(n->kids[0].get(), x, y, z, k);
                T b = evalExpr<T>(n->kids[1].get(), x, y, z, k);
                if (f == "pow") return mc_pow(a, b);
                if (f == "min") return mc_min(a, b);
                if (f == "max") return mc_max(a, b);
                throw ExprError("unknown function '" + f + "' (2 arguments)");
            }
            throw ExprError("function '" + f + "' called with unsupported argument count");
        }
    }
    throw ExprError("internal: unhandled node type");
}
