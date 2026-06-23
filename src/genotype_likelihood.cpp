/**
 * @file genotype_likelihood.cpp
 * @brief Genotype likelihood extraction and quality control implementation
 */

#include "genotype_likelihood.h"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace fastlckin {

std::vector<GenotypeLikelihood> extract_genotype_likelihoods(
    const ngslib::VCFRecord& rec,
    const ngslib::VCFHeader& hdr,
    int n_samples,
    int gq_min,
    const std::string& pl_field)
{
    std::vector<GenotypeLikelihood> result(n_samples);

    // Try the specified Phred-scaled field first (default: PL)
    std::vector<int32_t> pl_values;
    int pl_per_sample = rec.get_format_int(hdr, pl_field.c_str(), pl_values);

    // Try GL (float log10-scale) as fallback
    std::vector<float> gl_values;
    int gl_per_sample = 0;
    if (pl_per_sample < 3) {
        gl_per_sample = rec.get_format_float(hdr, "GL", gl_values);
    }

    // Try GQ
    std::vector<int32_t> gq_values;
    int gq_per_sample = rec.get_format_int(hdr, "GQ", gq_values);

    bool has_pl = (pl_per_sample >= 3);
    bool has_gl = (gl_per_sample >= 3);
    bool has_gq = (gq_per_sample >= 1);

    if (!has_pl && !has_gl) {
        // No likelihood data available: mask all
        for (int i = 0; i < n_samples; ++i) {
            result[i].masked = true;
        }
        return result;
    }

    for (int i = 0; i < n_samples; ++i) {
        auto& gl = result[i];

        if (has_pl) {
            int base = i * pl_per_sample;
            int32_t pl0 = pl_values[base + 0];
            int32_t pl1 = pl_values[base + 1];
            int32_t pl2 = pl_values[base + 2];

            // Check for missing values (bcf_int32_missing)
            if (pl0 == ngslib::VCFRecord::INT_MISSING ||
                pl1 == ngslib::VCFRecord::INT_MISSING ||
                pl2 == ngslib::VCFRecord::INT_MISSING) {
                gl.masked = true;
                continue;
            }

            // PL → linear likelihood: 10^(-PL/10)
            gl.gl[0] = std::pow(10.0, -static_cast<double>(pl0) / 10.0);
            gl.gl[1] = std::pow(10.0, -static_cast<double>(pl1) / 10.0);
            gl.gl[2] = std::pow(10.0, -static_cast<double>(pl2) / 10.0);

            // Derive GQ from PL if no GQ field: second smallest PL - smallest PL
            if (has_gq) {
                gl.gq = (gq_values[i] == ngslib::VCFRecord::INT_MISSING) ? 0 : gq_values[i];
            } else {
                int sorted_pl[3] = {pl0, pl1, pl2};
                std::sort(sorted_pl, sorted_pl + 3);
                gl.gq = sorted_pl[1] - sorted_pl[0];
            }
        } else {
            // GL → linear likelihood: 10^GL (GL is log10-scale)
            int base = i * gl_per_sample;
            float gl0 = gl_values[base + 0];
            float gl1 = gl_values[base + 1];
            float gl2 = gl_values[base + 2];

            // Check for missing float values
            if (std::isnan(gl0) || std::isnan(gl1) || std::isnan(gl2)) {
                gl.masked = true;
                continue;
            }

            gl.gl[0] = std::pow(10.0, static_cast<double>(gl0));
            gl.gl[1] = std::pow(10.0, static_cast<double>(gl1));
            gl.gl[2] = std::pow(10.0, static_cast<double>(gl2));

            // Derive GQ from GL: convert to PL-like, then compute
            double pl0 = -10.0 * gl0;
            double pl1 = -10.0 * gl1;
            double pl2 = -10.0 * gl2;
            double sorted_pl[3] = {pl0, pl1, pl2};
            std::sort(sorted_pl, sorted_pl + 3);
            gl.gq = static_cast<int>(sorted_pl[1] - sorted_pl[0]);

            if (has_gq) {
                gl.gq = (gq_values[i] == ngslib::VCFRecord::INT_MISSING) ? 0 : gq_values[i];
            }
        }

        // GQ quality filter
        if (gl.gq < gq_min) {
            gl.masked = true;
        }
    }

    return result;
}

}  // namespace fastlckin
