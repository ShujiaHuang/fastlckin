/**
 * @file optimizer.cpp
 * @brief Nelder-Mead simplex optimizer implementation
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include "optimizer.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace fastlckin {
namespace {

// NM standard parameters
constexpr double ALPHA   = 1.0;   // reflection
constexpr double GAMMA   = 2.0;   // expansion
constexpr double RHO     = 0.5;   // contraction
constexpr double SIGMA   = 0.5;   // shrink

/// Build initial simplex: n+1 vertices for n-dimensional problem.
/// Boundary-aware: when x0 is near the upper bound of the IBD simplex
/// (k0+k1+k2=1, i.e. x[0]+x[1] <= 1), steps are scaled to remain valid.
std::vector<std::vector<double>> make_simplex(const std::vector<double>& x0) {
    int n = static_cast<int>(x0.size());
    std::vector<std::vector<double>> simplex(n + 1, x0);

    // Sum of coordinates (for IBD: k1+k2; valid when <= 1)
    double x_sum = 0.0;
    for (int j = 0; j < n; ++j) x_sum += x0[j];

    for (int i = 0; i < n; ++i) {
        if (std::abs(x0[i]) > 1e-8) {
            double step = std::abs(x0[i]) * 0.5;
            // Don't step past the simplex boundary (sum <= 1)
            double remaining = std::max(0.01, 1.0 - x_sum);
            step = std::min(step, remaining);
            simplex[i + 1][i] += step;
        } else {
            double step = 0.5;
            double remaining = std::max(0.01, 1.0 - x_sum);
            simplex[i + 1][i] = std::min(step, remaining);
        }
    }
    return simplex;
}

/// Euclidean distance between two points
double point_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

}  // anonymous namespace

NelderMeadResult nelder_mead(
    std::function<double(const std::vector<double>&)> func,
    const std::vector<double>& x0,
    const NelderMeadConfig& config)
{
    int n = static_cast<int>(x0.size());
    NelderMeadResult result;
    result.x = x0;
    result.fval = func(x0);
    result.fun_evals = 1;

    // Build initial simplex and evaluate
    auto simplex = make_simplex(x0);
    std::vector<double> fvals(n + 1);
    for (int i = 0; i <= n; ++i) {
        fvals[i] = func(simplex[i]);
        ++result.fun_evals;
    }

    for (int iter = 0; iter < config.max_iter; ++iter) {
        result.iterations = iter + 1;

        // Sort vertices by function value
        std::vector<int> order(n + 1);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return fvals[a] < fvals[b]; });

        std::vector<std::vector<double>> sorted_simplex(n + 1);
        std::vector<double> sorted_fvals(n + 1);
        for (int i = 0; i <= n; ++i) {
            sorted_simplex[i] = simplex[order[i]];
            sorted_fvals[i] = fvals[order[i]];
        }
        simplex = sorted_simplex;
        fvals = sorted_fvals;

        // Check convergence
        double f_range = fvals[n] - fvals[0];
        double max_dist = 0.0;
        for (int i = 1; i <= n; ++i) {
            max_dist = std::max(max_dist, point_distance(simplex[0], simplex[i]));
        }

        if (f_range < config.ftol && max_dist < config.xtol) {
            result.converged = true;
            result.x = simplex[0];
            result.fval = fvals[0];
            return result;
        }

        if (result.fun_evals >= config.max_fun_evals) {
            result.x = simplex[0];
            result.fval = fvals[0];
            return result;
        }

        // Centroid of all points except the worst
        std::vector<double> centroid(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                centroid[j] += simplex[i][j];
            }
        }
        for (int j = 0; j < n; ++j) {
            centroid[j] /= n;
        }

        // Reflection
        std::vector<double> xr(n);
        for (int j = 0; j < n; ++j) {
            xr[j] = centroid[j] + ALPHA * (centroid[j] - simplex[n][j]);
        }
        double fr = func(xr);
        ++result.fun_evals;

        if (fr < fvals[0]) {
            // Expansion
            std::vector<double> xe(n);
            for (int j = 0; j < n; ++j) {
                xe[j] = centroid[j] + GAMMA * (xr[j] - centroid[j]);
            }
            double fe = func(xe);
            ++result.fun_evals;

            if (fe < fr) {
                simplex[n] = xe;
                fvals[n] = fe;
            } else {
                simplex[n] = xr;
                fvals[n] = fr;
            }
        } else if (fr < fvals[n - 1]) {
            // Accept reflection
            simplex[n] = xr;
            fvals[n] = fr;
        } else {
            // Contraction
            std::vector<double> xc(n);
            if (fr < fvals[n]) {
                // Outside contraction
                for (int j = 0; j < n; ++j) {
                    xc[j] = centroid[j] + RHO * (xr[j] - centroid[j]);
                }
            } else {
                // Inside contraction
                for (int j = 0; j < n; ++j) {
                    xc[j] = centroid[j] + RHO * (simplex[n][j] - centroid[j]);
                }
            }
            double fc = func(xc);
            ++result.fun_evals;

            if (fc < std::min(fr, fvals[n])) {
                simplex[n] = xc;
                fvals[n] = fc;
            } else {
                // Shrink
                for (int i = 1; i <= n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        simplex[i][j] = simplex[0][j]
                                        + SIGMA * (simplex[i][j] - simplex[0][j]);
                    }
                    fvals[i] = func(simplex[i]);
                    ++result.fun_evals;
                }
            }
        }
    }

    // Final sort to get best
    int best = 0;
    for (int i = 1; i <= n; ++i) {
        if (fvals[i] < fvals[best]) best = i;
    }
    result.x = simplex[best];
    result.fval = fvals[best];
    return result;
}

}  // namespace fastlckin
