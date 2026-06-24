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
    std::vector<double> standard_errors;  ///< Standard errors for each parameter (empty if not computed)
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

/// Compute standard errors via numerical Hessian at the MLE point.
///
/// Uses finite difference approximation to compute the Hessian matrix,
/// then inverts it to get the covariance matrix. Standard errors are
/// the square roots of the diagonal elements.
///
/// @param func   Objective function (negative log-likelihood)
/// @param x_opt  Optimal parameter values
/// @param h      Step size for finite difference (default: 1e-5)
/// @return Vector of standard errors, one per parameter
std::vector<double> compute_standard_errors(
    std::function<double(const std::vector<double>&)> func,
    const std::vector<double>& x_opt,
    double h = 1e-6
);

}  // namespace fastlckin

#endif
