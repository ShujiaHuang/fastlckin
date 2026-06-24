#!/usr/bin/env python3
"""
validate_ground_truth.py
========================
Validate fastlckin kinship estimation against known-relationship synthetic data.

Workflow:
  1. Generate synthetic PLINK data with known relationships (if not present)
  2. Run fastlckin on the data
  3. Compare estimated IBD coefficients (k0, k1, k2, PI_HAT) with ground truth
  4. Report accuracy per relationship type and overall PASS/FAIL

Usage:
  cd <project_root>
  python3 tests/validate_ground_truth.py [--fastlckin PATH] [--data-dir DIR]
"""

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────

# Tolerances for each relationship type (|estimated - expected| <= tol)
# These account for sampling variance with ~500 SNPs.
TOLERANCES = {
    "Duplicate/MZ":     {"k0": 0.05, "k1": 0.05, "k2": 0.05, "pi_hat": 0.05},
    "Parent-Offspring": {"k0": 0.15, "k1": 0.20, "k2": 0.15, "pi_hat": 0.15},
    "Full-Sibling":     {"k0": 0.20, "k1": 0.25, "k2": 0.20, "pi_hat": 0.20},
    "Half-Sibling":     {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Avuncular":        {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Grandparent":      {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Unrelated":        {"k0": 0.15, "k1": 0.15, "k2": 0.05, "pi_hat": 0.10},
}

# Color codes for terminal output
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"

# ── Helpers ───────────────────────────────────────────────────────────


def normalize_pair(ind1: str, ind2: str) -> tuple:
    """Return a canonical (sorted) pair key."""
    return tuple(sorted([ind1, ind2]))


def load_ground_truth(path: str) -> dict:
    """Load ground_truth.tsv into a dict keyed by normalized pair."""
    gt = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = normalize_pair(row["Ind1"], row["Ind2"])
            gt[key] = {
                "relationship": row["Relationship"],
                "k0": float(row["Expected_k0"]),
                "k1": float(row["Expected_k1"]),
                "k2": float(row["Expected_k2"]),
                "pi_hat": float(row["Expected_PI_HAT"]),
            }
    return gt


def load_fastlckin_output(path: str) -> dict:
    """Load fastlckin output TSV into a dict keyed by normalized pair."""
    results = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Ind1"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            ind1, ind2 = parts[0], parts[1]
            key = normalize_pair(ind1, ind2)
            results[key] = {
                "k0": float(parts[2]),
                "k1": float(parts[3]),
                "k2": float(parts[4]),
                "pi_hat": float(parts[5]),
                "n_snps": int(parts[6]),
                "relationship": parts[7] if len(parts) > 7 else "",
            }
    return results


def check_tolerance(est: float, exp: float, tol: float) -> bool:
    """Check if |est - exp| <= tol."""
    return abs(est - exp) <= tol


def fmt_val(est: float, exp: float, tol: float, passed: bool) -> str:
    """Format a value with color coding."""
    color = GREEN if passed else RED
    sign = "✓" if passed else "✗"
    return f"{color}{sign} {est:.4f} (exp={exp:.3f}±{tol}){RESET}"


# ── Main validation logic ────────────────────────────────────────────


def run_validation(fastlckin_path: str, data_dir: str, prefix: str,
                   skip_generation: bool, verbose: bool) -> int:
    """Run the full validation pipeline. Returns exit code (0=pass, 1=fail)."""

    plink_prefix = os.path.join(data_dir, prefix)
    gt_path = os.path.join(data_dir, "ground_truth.tsv")
    output_path = os.path.join(data_dir, "kinship_validation_output.tsv")

    # ── Step 1: Generate data ─────────────────────────────────────────
    if not skip_generation:
        gen_bin = None
        # Look for the generator in build/tests/
        for candidate in [
            "build/tests/generate_known_relationships",
            "bin/generate_known_relationships",
        ]:
            if os.path.isfile(candidate):
                gen_bin = candidate
                break

        if gen_bin is None:
            print(f"{RED}ERROR: generate_known_relationships binary not found.{RESET}")
            print("  Build it first: cd build && make generate_known_relationships")
            return 1

        print(f"{CYAN}[1/4] Generating synthetic data with known relationships...{RESET}")
        result = subprocess.run(
            [gen_bin], capture_output=True, text=True, cwd=os.getcwd()
        )
        if result.returncode != 0:
            print(f"{RED}ERROR: Data generation failed:{RESET}")
            print(result.stderr)
            return 1
        if verbose:
            print(result.stderr)
    else:
        print(f"{CYAN}[1/4] Skipping data generation (using existing data).{RESET}")

    # Verify files exist
    for ext in [".bed", ".bim", ".fam"]:
        fpath = plink_prefix + ext
        if not os.path.isfile(fpath):
            print(f"{RED}ERROR: Missing {fpath}{RESET}")
            return 1
    if not os.path.isfile(gt_path):
        print(f"{RED}ERROR: Missing {gt_path}{RESET}")
        return 1

    # ── Step 2: Run fastlckin ─────────────────────────────────────────
    print(f"{CYAN}[2/4] Running fastlckin relatedness estimation...{RESET}")
    cmd = [
        fastlckin_path, "relatedness",
        "-p", plink_prefix,
        "--no-ld-prune",
        "--classify",
        "-o", output_path,
    ]
    if verbose:
        cmd.append("--verbose")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"{RED}ERROR: fastlckin failed:{RESET}")
        print(result.stderr)
        return 1
    print(f"  {result.stderr.strip()}")

    # ── Step 3: Load and compare ──────────────────────────────────────
    print(f"{CYAN}[3/4] Comparing estimated vs. expected IBD coefficients...{RESET}")

    gt = load_ground_truth(gt_path)
    est = load_fastlckin_output(output_path)

    print(f"  Ground truth pairs : {len(gt)}")
    print(f"  Estimated pairs   : {len(est)}")
    print()

    # ── Step 4: Detailed comparison ───────────────────────────────────
    print(f"{CYAN}[4/4] Detailed results:{RESET}")
    print("=" * 110)
    print(f"{'Pair':<28} {'Type':<18} {'k0':>28} {'k1':>28} {'k2':>28}")
    print("-" * 110)

    n_pass = 0
    n_fail = 0
    n_missing = 0
    stats_by_type = defaultdict(lambda: {"pass": 0, "fail": 0})

    for key, expected in sorted(gt.items(), key=lambda x: x[1]["relationship"]):
        pair_str = f"{key[0]:<12} {key[1]:<12}"
        rel = expected["relationship"]
        tol = TOLERANCES.get(rel, TOLERANCES["Unrelated"])

        if key not in est:
            print(f"  {pair_str}  {rel:<18} {YELLOW}MISSING{RESET}")
            n_missing += 1
            stats_by_type[rel]["fail"] += 1
            continue

        e = est[key]

        checks = {
            "k0": check_tolerance(e["k0"], expected["k0"], tol["k0"]),
            "k1": check_tolerance(e["k1"], expected["k1"], tol["k1"]),
            "k2": check_tolerance(e["k2"], expected["k2"], tol["k2"]),
            "pi_hat": check_tolerance(e["pi_hat"], expected["pi_hat"], tol["pi_hat"]),
        }

        all_pass = all(checks.values())
        if all_pass:
            n_pass += 1
            stats_by_type[rel]["pass"] += 1
        else:
            n_fail += 1
            stats_by_type[rel]["fail"] += 1

        status = f"{GREEN}PASS{RESET}" if all_pass else f"{RED}FAIL{RESET}"
        k0_s = fmt_val(e["k0"], expected["k0"], tol["k0"], checks["k0"])
        k1_s = fmt_val(e["k1"], expected["k1"], tol["k1"], checks["k1"])
        k2_s = fmt_val(e["k2"], expected["k2"], tol["k2"], checks["k2"])

        print(f"  {pair_str}  {rel:<18} {k0_s:>40} {k1_s:>40} {k2_s:>40}")

        if not all_pass or verbose:
            pi_s = fmt_val(e["pi_hat"], expected["pi_hat"], tol["pi_hat"],
                           checks["pi_hat"])
            print(f"  {'':>30} PI_HAT: {pi_s}  [{status}]")

    print("=" * 110)

    # ── Summary by relationship type ──────────────────────────────────
    print(f"\n{BOLD}Summary by relationship type:{RESET}")
    print(f"  {'Type':<20} {'Pass':>6} {'Fail':>6} {'Total':>6} {'Rate':>8}")
    print(f"  {'-'*50}")
    for rel in sorted(stats_by_type.keys()):
        s = stats_by_type[rel]
        total = s["pass"] + s["fail"]
        rate = s["pass"] / total * 100 if total > 0 else 0
        color = GREEN if s["fail"] == 0 else RED
        print(f"  {rel:<20} {s['pass']:>6} {s['fail']:>6} {total:>6} "
              f"{color}{rate:>7.1f}%{RESET}")

    # ── Overall summary ───────────────────────────────────────────────
    total_checked = n_pass + n_fail
    overall_rate = n_pass / total_checked * 100 if total_checked > 0 else 0
    print(f"\n{BOLD}Overall:{RESET} {n_pass}/{total_checked} pairs passed "
          f"({overall_rate:.1f}%), {n_missing} missing")

    if n_fail == 0 and n_missing == 0:
        print(f"\n{GREEN}{BOLD}✓ ALL VALIDATIONS PASSED{RESET}")
        return 0
    else:
        print(f"\n{RED}{BOLD}✗ SOME VALIDATIONS FAILED{RESET}")
        print(f"  {n_fail} pair(s) exceeded tolerance, {n_missing} pair(s) missing.")
        print(f"\n  Note: Tolerances are set for ~500 SNPs with random allele")
        print(f"  frequencies. With few SNPs, estimates have inherent sampling")
        print(f"  variance. Review failed pairs and consider adjusting tolerances")
        print(f"  or increasing SNP count in the data generator.")
        return 1


# ── Entry point ───────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Validate fastlckin against known-relationship ground truth"
    )
    parser.add_argument(
        "--fastlckin", default="bin/fastlckin",
        help="Path to fastlckin binary (default: bin/fastlckin)"
    )
    parser.add_argument(
        "--data-dir", default="tests/data",
        help="Directory for test data (default: tests/data)"
    )
    parser.add_argument(
        "--prefix", default="kinship_test",
        help="PLINK file prefix (default: kinship_test)"
    )
    parser.add_argument(
        "--skip-generation", action="store_true",
        help="Skip data generation (use existing files)"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Show detailed output"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.fastlckin):
        print(f"{RED}ERROR: fastlckin binary not found at {args.fastlckin}{RESET}")
        print("  Build it first: cd build && make fastlckin")
        sys.exit(1)

    exit_code = run_validation(
        fastlckin_path=args.fastlckin,
        data_dir=args.data_dir,
        prefix=args.prefix,
        skip_generation=args.skip_generation,
        verbose=args.verbose,
    )
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
