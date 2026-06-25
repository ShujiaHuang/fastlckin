#ifndef _FASTLCKIN_OPTIMIZER_H_
#define _FASTLCKIN_OPTIMIZER_H_

/**
 * @file optimizer.h
 * @brief Nelder-Mead simplex optimizer
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include <vector>
#include <functional>

namespace fastlckin {

struct NelderMeadConfig {
    double xtol = 1e-4;         ///< Parameter convergence threshold (improved from 0.01 in v0.3.0)
    double ftol = 1e-4;         ///< Function value convergence threshold (improved from 0.01 in v0.3.0)
    int max_iter = 10000;       ///< Maximum iterations
    int max_fun_evals = 10000;  ///< Maximum function evaluations
};

struct NelderMeadResult {
    std::vector<double> x;      ///< Optimal point
    double fval = 0.0;          ///< Optimal function value
    int iterations = 0;         ///< Number of iterations performed
    int fun_evals = 0;          ///< Number of function evaluations
    bool converged = false;     ///< Whether optimization converged
};

/// Nelder-Mead simplex method (minimization)
/// @param func  Objective function f(x) -> double
/// @param x0    Initial point (2D: [k1, k2])
/// @param config Optimization parameters
/// @return Optimization result
NelderMeadResult nelder_mead(
    std::function<double(const std::vector<double>&)> func,
    const std::vector<double>& x0,
    const NelderMeadConfig& config = NelderMeadConfig()
);

}  // namespace fastlckin

#endif
