#!/bin/bash
# run_eval.sh — fastlckin evaluation pipeline
# Usage: bash tests/run_eval.sh [--skip-gen]
#
# ═══════════════════════════════════════════════════════════════════════
#  Evaluation Overview
# ═══════════════════════════════════════════════════════════════════════
#
#  This script runs fastlckin through multiple evaluation scenarios:
#
#  Step 1-3: Large-scale evaluation (5k SNPs, 4 modes)
#    - Tests IBD coefficient accuracy (k0, k1, k2, PI_HAT)
#    - Tests classification accuracy (threshold vs likelihood methods)
#    - Compares PLINK (hard genotypes) vs VCF (genotype likelihoods)
#
#  Step 4: Small-scale validation (20k SNPs, PLINK-only)
#    - High-SNP-count validation with --no-ld-prune
#    - Verifies IBD accuracy on a different dataset
#
#  Expected results if everything is correct:
#    - IBD accuracy: ~100% pairs within tolerance
#    - PI_HAT for unrelated pairs: ≈0 (typically <0.02)
#    - PI_HAT for parent-offspring: ≈0.5
#    - PI_HAT for full-siblings: ≈0.5
#    - PI_HAT for 2nd-degree (half-sib, avuncular, grandparent): ≈0.25
#    - PI_HAT for duplicate/MZ: ≈1.0
#    - Threshold classification: ~100% for coarse categories (UN/DUP/PO/FS/2nd)
#    - Likelihood classification: ~85% fine-grained (2nd-degree subtypes hard to distinguish)
#
# ═══════════════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

FASTLCKIN="$PROJECT_DIR/bin/fastlckin"
EVAL_DIR=/tmp/eval_data
EVALUATE="python3 $SCRIPT_DIR/evaluate.py"

# ── Step 1: Generate 5k evaluation dataset ────────────────────────────
# Generates synthetic data with known relationships:
#   - 60 founders, 5000 SNPs
#   - Relationship types: Duplicate/MZ, Parent-Offspring, Full-Sibling,
#     Half-Sibling, Avuncular, Grandparent, Unrelated
#   - Also creates VCF with simulated sequencing depth ~8x
#
# Output:
#   - $EVAL_DIR/eval_5k.bed/bim/fam (PLINK files)
#   - $EVAL_DIR/eval_5k.vcf.gz (VCF with PL field)
#   - $EVAL_DIR/ground_truth.tsv (known relationships)
if [[ "${1:-}" != "--skip-gen" ]]; then
    echo "=== Generating 5k evaluation data ==="
    "$PROJECT_DIR/build/tests/generate_known_relationships" \
        --founders 60 --snps 5000 --dup 3 --po 3 --fs 5 --hs 5 \
        --af-min 0.1 --af-max 0.5 --seed 42 \
        --prefix eval_5k --output-dir "$EVAL_DIR"
    python3 "$SCRIPT_DIR/plink_to_vcf_pl.py" "$EVAL_DIR/eval_5k" "$EVAL_DIR/eval_5k.vcf" \
        --depth 8 --seed 42
    bgzip -f "$EVAL_DIR/eval_5k.vcf"
    tabix -f -p vcf "$EVAL_DIR/eval_5k.vcf.gz"
else
    echo "=== Skipping data generation (--skip-gen) ==="
fi

# ── Step 2: Run fastlckin in all modes ───────────────────────────────
# Tests 4 combinations:
#   1. PLINK + Threshold: Hard genotypes, rule-based classification
#   2. PLINK + Likelihood: Hard genotypes, model-selection classification
#   3. VCF + Threshold: Genotype likelihoods (uncertain), rule-based
#   4. VCF + Likelihood: Genotype likelihoods (uncertain), model-selection
#
# Expected: All modes should produce similar IBD estimates
#           VCF modes may be slightly less accurate due to GL uncertainty
echo "=== PLINK Threshold mode ==="
$FASTLCKIN relatedness --threads 10 --classify --output-likelihoods \
    -p "$EVAL_DIR/eval_5k" -o "$EVAL_DIR/threshold.tsv"

echo "=== PLINK Likelihood mode ==="
$FASTLCKIN relatedness --threads 10 --classify=likelihood --output-likelihoods \
    -p "$EVAL_DIR/eval_5k" -o "$EVAL_DIR/likelihood.tsv"

echo "=== VCF Threshold mode ==="
$FASTLCKIN relatedness --threads 10 --classify --output-likelihoods --threads 16 \
    -v "$EVAL_DIR/eval_5k.vcf.gz" -o "$EVAL_DIR/vcf_threshold.tsv"

echo "=== VCF Likelihood mode ==="
$FASTLCKIN relatedness --threads 10 --classify=likelihood --output-likelihoods --threads 16 \
    -v "$EVAL_DIR/eval_5k.vcf.gz" -o "$EVAL_DIR/vcf_likelihood.tsv"

# ── Step 3: Evaluate ─────────────────────────────────────────────────
# Compares fastlckin output against ground truth:
#   - IBD accuracy: k0, k1, k2, PI_HAT within tolerance
#   - PI_HAT statistics: mean/std per relationship type
#   - Classification: confusion matrices for threshold/likelihood methods
#
# Expected results (if correct):
#   - IBD: ~94/94 pairs pass (100%)
#   - Unrelated PI_HAT: mean < 0.02
#   - Threshold coarse classification: ~94/94 (100%) for PLINK-TH, ~93/94 for VCF-TH
#   - Likelihood fine classification: ~80/94 (85%), errors in 2nd-degree subtypes
echo ""
echo "=========================================="
echo "  Evaluation Report (5k dataset)"
echo "=========================================="

$EVALUATE \
    --label PLINK-TH --label PLINK-LK --label VCF-TH --label VCF-LK \
    "$EVAL_DIR/ground_truth.tsv" \
    "$EVAL_DIR/threshold.tsv" \
    "$EVAL_DIR/likelihood.tsv" \
    "$EVAL_DIR/vcf_threshold.tsv" \
    "$EVAL_DIR/vcf_likelihood.tsv"

# ── Step 4: Small-scale validation (tests/data) ──────────────────────
# Uses kit_20k dataset:
#   - 20000 SNPs (higher than 5k eval)
#   - PLINK-only mode with --no-ld-prune (all SNPs used)
#   - Tests IBD accuracy on a different sample composition
#
# Expected results (if correct):
#   - IBD: ~23/23 pairs pass (100%) — more SNPs = better accuracy
#   - Unrelated PI_HAT: ≈0
#   - Threshold classification: high accuracy for major categories
KIT20K_DIR="$SCRIPT_DIR/data"
if [[ ! -f "$KIT20K_DIR/kit_20k.bed" ]]; then
    echo ""
    echo "=========================================="
    echo "  Generating kit_20k test data"
    echo "=========================================="
    "$PROJECT_DIR/build/tests/generate_known_relationships" \
        --snps 20000 --seed 42 --prefix kit_20k --output-dir "$KIT20K_DIR"
    # Rename ground_truth.tsv to avoid conflict with 5k eval data
    mv "$KIT20K_DIR/ground_truth.tsv" "$KIT20K_DIR/kit_20k_ground_truth.tsv"
fi

echo ""
echo "=========================================="
echo "  Small-scale Validation (kit_20k, 20k SNPs)"
echo "=========================================="

$FASTLCKIN relatedness --threads 10 -p "$KIT20K_DIR/kit_20k" --no-ld-prune --classify \
    -o "$KIT20K_DIR/kit_20k_output.tsv"

$EVALUATE "$KIT20K_DIR/kit_20k_ground_truth.tsv" "$KIT20K_DIR/kit_20k_output.tsv"

echo ""
echo "=========================================="
echo "  Evaluation Complete"
echo "=========================================="
echo "If all tests passed, you should see:"
echo "  - IBD accuracy: ~100% for both datasets"
echo "  - Unrelated PI_HAT: ≈0 (<0.02)"
echo "  - Classification: high accuracy for major relationship types"
