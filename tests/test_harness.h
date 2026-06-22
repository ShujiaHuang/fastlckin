#ifndef _FASTLCKIN_TEST_HARNESS_H_
#define _FASTLCKIN_TEST_HARNESS_H_

/**
 * @file test_harness.h
 * @brief Minimal unit-test framework (no external dependencies)
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <sstream>

namespace test {

struct TestCase {
    std::string name;
    std::function<void()> func;
};

static std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

static int g_pass = 0;
static int g_fail = 0;

struct Registrar {
    Registrar(const std::string& name, std::function<void()> f) {
        registry().push_back({name, std::move(f)});
    }
};

inline void check_impl(bool cond, const char* expr, const char* file, int line) {
    if (cond) {
        ++g_pass;
    } else {
        ++g_fail;
        std::cerr << "  FAIL: " << expr << "  (" << file << ":" << line << ")\n";
    }
}

inline void check_near_impl(double a, double b, double tol, const char* expr, const char* file, int line) {
    if (std::abs(a - b) <= tol) {
        ++g_pass;
    } else {
        ++g_fail;
        std::cerr << "  FAIL: " << expr << " (got " << a << ", expected " << b
                  << ", tol=" << tol << ")  (" << file << ":" << line << ")\n";
    }
}

#define TEST_CASE(name) \
    static void test_##name(); \
    static test::Registrar reg_##name(#name, test_##name); \
    static void test_##name()

#define CHECK(expr) test::check_impl((expr), #expr, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol) test::check_near_impl((a), (b), (tol), #a " ~= " #b, __FILE__, __LINE__)

inline int run_all() {
    g_pass = 0;
    g_fail = 0;
    for (auto& tc : registry()) {
        std::cerr << "[TEST] " << tc.name << "\n";
        try {
            tc.func();
        } catch (const std::exception& e) {
            ++g_fail;
            std::cerr << "  EXCEPTION: " << e.what() << "\n";
        }
    }
    std::cerr << "\n===== Results: " << g_pass << " passed, " << g_fail << " failed =====\n";
    return (g_fail > 0) ? 1 : 0;
}

}  // namespace test

#define RUN_ALL_TESTS() test::run_all()

#endif
