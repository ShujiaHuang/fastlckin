#!/usr/bin/env python3
"""
evaluate.py — fastlckin evaluation tool
========================================

Compare fastlckin output against ground truth and report:
  - IBD coefficient accuracy (k0, k1, k2, PI_HAT) per relationship type
  - Classification confusion matrices (if classification columns present)
  - PI_HAT statistics per relationship type

Usage:
  python3 tests/evaluate.py <ground_truth> <result> [result2 ...]
  python3 tests/evaluate.py --label PLINK --label VCF gt.tsv plink_out.tsv vcf_out.tsv

Arguments:
  ground_truth   Path to ground_truth.tsv (with Ind1, Ind2, Relationship, Expected_*)
  result         One or more fastlckin output TSV files

Examples:
  # Single result
  python3 tests/evaluate.py tests/data/ground_truth.tsv tests/data/kinship_validation_output.tsv

  # Multiple results for comparison
  python3 tests/evaluate.py /tmp/eval_data/ground_truth.tsv \\
      /tmp/eval_data/threshold.tsv /tmp/eval_data/vcf_threshold.tsv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

# ═══════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════

RED    = "\033[91m"
GREEN  = "\033[92m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

EXPECTED_PIHAT = {
    'Duplicate/MZ': 1.0, 'Parent-Offspring': 0.5, 'Full-Sibling': 0.5,
    'Half-Sibling': 0.25, 'Avuncular': 0.25, 'Grandparent': 0.25,
    'Great-Grandparent': 0.125, 'Unrelated': 0.0,
}

COARSE_MAP = {
    'Duplicate/MZ': 'DUP', 'Parent-Offspring': 'PO', 'Full-Sibling': 'FS',
    'Half-Sibling': '2nd', 'Avuncular': '2nd', 'Grandparent': '2nd',
    'Great-Grandparent': '3rd', 'Unrelated': 'UN',
}
THRESHOLD_LABELS = {
    'Duplicate/MZ': 'DUP', 'Parent-Offspring': 'PO', 'Full-Sibling': 'FS',
    'Second-degree': '2nd', 'Third-degree': '3rd', 'Unrelated': 'UN',
}
FINE_MAP = {
    # Ground truth relationships map to IBD equivalence classes
    'Duplicate/MZ': 'Duplicate/MZ', 'Parent-Offspring': 'Parent-Offspring',
    'Full-Sibling': 'Full-Sibling', 'Half-Sibling': '2nd-Degree',
    'Avuncular': '2nd-Degree', 'Grandparent': '2nd-Degree',
    'First-Cousin': '3rd-Degree', 'Great-Grandparent': '3rd-Degree',
    'Unrelated': 'Unrelated',
    # Predictions already use IBD equivalence class names
    '2nd-Degree': '2nd-Degree', '3rd-Degree': '3rd-Degree',
}

TOLERANCES = {
    "Duplicate/MZ":     {"k0": 0.05, "k1": 0.05, "k2": 0.05, "pi_hat": 0.05},
    "Parent-Offspring": {"k0": 0.15, "k1": 0.20, "k2": 0.15, "pi_hat": 0.15},
    "Full-Sibling":     {"k0": 0.20, "k1": 0.25, "k2": 0.20, "pi_hat": 0.20},
    "2nd-Degree":       {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "3rd-Degree":       {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Unrelated":        {"k0": 0.15, "k1": 0.15, "k2": 0.05, "pi_hat": 0.10},
    # Legacy names for backward compatibility
    "Half-Sibling":     {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Avuncular":        {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Grandparent":      {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "First-Cousin":     {"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
    "Great-Grandparent":{"k0": 0.25, "k1": 0.25, "k2": 0.15, "pi_hat": 0.15},
}


# ═══════════════════════════════════════════════════════════════════════
#  Data loading
# ═══════════════════════════════════════════════════════════════════════

def normalize_pair(ind1: str, ind2: str) -> tuple:
    return tuple(sorted([ind1, ind2]))


def load_ground_truth(path: str) -> dict:
    """Load ground_truth.tsv → {(ind1,ind2): {'relationship': str, 'k0': float, ...}}"""
    gt = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = normalize_pair(row["Ind1"], row["Ind2"])
            entry = {"relationship": row["Relationship"]}
            for col, out in [("Expected_k0", "k0"), ("Expected_k1", "k1"),
                             ("Expected_k2", "k2"), ("Expected_PI_HAT", "pi_hat")]:
                if col in row:
                    entry[out] = float(row[col])
            gt[key] = entry
    return gt


def load_results(path: str) -> dict:
    """Load fastlckin output TSV → {(ind1,ind2): {col: value, ...}}"""
    results = {}
    with open(path) as f:
        header = None
        for line in f:
            if not line.startswith("#"):
                header = line.strip().split("\t")
                break
        if not header:
            return results
        col_idx = {h: i for i, h in enumerate(header)}
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            key = normalize_pair(cols[0], cols[1])
            results[key] = {h: cols[col_idx[h]] for h in col_idx if col_idx[h] < len(cols)}
    return results


# ═══════════════════════════════════════════════════════════════════════
#  Analysis: IBD accuracy
# ═══════════════════════════════════════════════════════════════════════

def _check(est, exp, tol):
    return abs(est - exp) <= tol


def _fmt(est, exp, tol, ok):
    c = GREEN if ok else RED
    s = "✓" if ok else "✗"
    return f"{c}{s} {est:.4f} (exp={exp:.3f}±{tol}){RESET}"


def report_ibd_accuracy(gt, results, label, out=sys.stdout):
    """Report per-pair IBD coefficient accuracy with tolerance checks."""
    col_map = {"k0": "k0", "k1": "k1", "k2": "k2", "pi_hat": "PI_HAT"}

    print(f"\n{'='*110}", file=out)
    print(f"  IBD Accuracy — {label}", file=out)
    print(f"{'='*110}", file=out)
    print(f"{'Pair':<28} {'Type':<20} {'k0':>30} {'k1':>30} {'k2':>30}", file=out)
    print(f"{'-'*110}", file=out)

    n_pass = n_fail = n_missing = 0
    stats = defaultdict(lambda: {"pass": 0, "fail": 0})

    for key, exp in sorted(gt.items(), key=lambda x: x[1]["relationship"]):
        pair_str = f"{key[0]:<12} {key[1]:<12}"
        rel = exp["relationship"]
        tol = TOLERANCES.get(rel, TOLERANCES["Unrelated"])

        if key not in results:
            print(f"  {pair_str}  {rel:<20} {YELLOW}MISSING{RESET}", file=out)
            n_missing += 1
            stats[rel]["fail"] += 1
            continue

        r = results[key]
        checks = {}
        for m in ["k0", "k1", "k2", "pi_hat"]:
            col = col_map[m]
            if m in exp and col in r:
                checks[m] = _check(float(r[col]), exp[m], tol[m])
            else:
                checks[m] = True

        ok = all(checks.values())
        if ok:
            n_pass += 1; stats[rel]["pass"] += 1
        else:
            n_fail += 1; stats[rel]["fail"] += 1

        k0_s = _fmt(float(r.get("k0", 0)), exp.get("k0", 0), tol["k0"], checks["k0"])
        k1_s = _fmt(float(r.get("k1", 0)), exp.get("k1", 0), tol["k1"], checks["k1"])
        k2_s = _fmt(float(r.get("k2", 0)), exp.get("k2", 0), tol["k2"], checks["k2"])
        print(f"  {pair_str}  {rel:<20} {k0_s:>32} {k1_s:>32} {k2_s:>32}", file=out)

    # Summary
    print(f"\n  {'Type':<20} {'Pass':>6} {'Fail':>6} {'Total':>6} {'Rate':>8}", file=out)
    print(f"  {'-'*50}", file=out)
    for rel in sorted(stats):
        s = stats[rel]
        t = s["pass"] + s["fail"]
        rate = s["pass"] / t * 100 if t else 0
        c = GREEN if s["fail"] == 0 else RED
        print(f"  {rel:<20} {s['pass']:>6} {s['fail']:>6} {t:>6} {c}{rate:>7.1f}%{RESET}", file=out)

    total = n_pass + n_fail
    rate = n_pass / total * 100 if total else 0
    print(f"\n  Overall: {n_pass}/{total} passed ({rate:.1f}%), {n_missing} missing", file=out)
    return n_fail, n_missing


# ═══════════════════════════════════════════════════════════════════════
#  Analysis: PI_HAT statistics
# ═══════════════════════════════════════════════════════════════════════

def report_pihat_stats(gt, results, label, out=sys.stdout):
    """Report PI_HAT mean/std per relationship type."""
    stats = defaultdict(list)
    for key in gt:
        if key in results and "PI_HAT" in results[key]:
            rel = gt[key]["relationship"]
            stats[rel].append(float(results[key]["PI_HAT"]))

    print(f"\n{'='*70}", file=out)
    print(f"  PI_HAT Statistics — {label}", file=out)
    print(f"{'='*70}", file=out)
    print(f"{'Relationship':<22}{'N':>5}{'Mean':>10}{'Std':>10}{'Expected':>10}", file=out)
    print(f"{'-'*57}", file=out)
    for rel in sorted(stats):
        vals = stats[rel]
        n = len(vals)
        mean = sum(vals) / n
        std = (sum((v - mean) ** 2 for v in vals) / n) ** 0.5
        print(f"{rel:<22}{n:>5}{mean:>10.4f}{std:>10.4f}{EXPECTED_PIHAT.get(rel, '?'):>10}", file=out)


# ═══════════════════════════════════════════════════════════════════════
#  Analysis: Classification
# ═══════════════════════════════════════════════════════════════════════

def _build_confusion(gt, pred, gt_map, pred_map, pred_col):
    matrix = defaultdict(lambda: defaultdict(int))
    correct = total = 0
    cats = set()
    for key in gt:
        if key not in pred:
            continue
        gt_rel = gt[key]["relationship"]
        gt_cat = gt_map.get(gt_rel, "Other")
        pred_val = pred[key].get(pred_col, "Unrelated")
        pred_cat = pred_map.get(pred_val, pred_val)
        matrix[gt_cat][pred_cat] += 1
        total += 1
        if gt_cat == pred_cat:
            correct += 1
        cats.update([gt_cat, pred_cat])
    return matrix, correct, total, sorted(cats)


def _print_confusion(matrix, cats, title, out=sys.stdout):
    w = max(max(len(c) for c in cats), 6) + 2
    print(f"\n{'='*70}", file=out)
    print(f"  {title}", file=out)
    print(f"{'='*70}", file=out)
    print(f"{'GT \\ Pred':<{w}}", end="", file=out)
    for c in cats:
        print(f"{c:>{w}}", end="", file=out)
    print(f"{'Total':>{w}}{'Recall':>{w}}", file=out)
    for gc in cats:
        print(f"{gc:<{w}}", end="", file=out)
        rt = rc = 0
        for pc in cats:
            v = matrix[gc][pc]; rt += v
            if gc == pc: rc = v
            print(f"{v:>{w}}", end="", file=out)
        print(f"{rt:>{w}}{rc/rt*100 if rt else 0:>{w-1}.1f}%", file=out)
    print(f"{'Total':<{w}}", end="", file=out)
    gt_total = 0
    for pc in cats:
        ct = sum(matrix[gc][pc] for gc in cats); gt_total += ct
        print(f"{ct:>{w}}", end="", file=out)
    print(f"{gt_total:>{w}}", file=out)


def report_classification(gt, results, label, out=sys.stdout):
    """Report classification confusion matrices and misclassifications."""
    has_threshold = any("Relationship" in results.get(k, {}) for k in list(results.keys())[:5])
    has_likelihood = any("MS_Relationship" in results.get(k, {}) for k in list(results.keys())[:5])

    summary = []

    if has_threshold:
        pred_coarse = {v: v for v in THRESHOLD_LABELS.values()}
        m, c, t, cats = _build_confusion(
            gt, results, COARSE_MAP, {**THRESHOLD_LABELS, **pred_coarse}, "Relationship")
        _print_confusion(m, cats, f"Threshold (coarse) — {label} — {c}/{t} = {c/t*100:.1f}%", out)
        summary.append(("Threshold (coarse)", c, t))

    if has_likelihood:
        m, c, t, cats = _build_confusion(gt, results, FINE_MAP, FINE_MAP, "MS_Relationship")
        _print_confusion(m, cats, f"Likelihood (fine) — {label} — {c}/{t} = {c/t*100:.1f}%", out)

        # Misclassifications
        errs = []
        for key in sorted(gt):
            if key not in results: continue
            gt_rel = gt[key]["relationship"]
            pred_rel = results[key].get("MS_Relationship", "Unrelated")
            if FINE_MAP.get(gt_rel) != FINE_MAP.get(pred_rel, pred_rel):
                errs.append((key[0], key[1], gt_rel, pred_rel,
                             float(results[key].get("k0", 0)),
                             float(results[key].get("PI_HAT", 0))))
        if errs:
            print(f"\n  Misclassifications:", file=out)
            print(f"  {'Ind1':<15}{'Ind2':<15}{'GT':<22}{'Pred':<22}{'k0':>8}{'PI_HAT':>8}", file=out)
            print(f"  {'-'*90}", file=out)
            for i1, i2, g, p, k0, ph in errs:
                print(f"  {i1:<15}{i2:<15}{g:<22}{p:<22}{k0:>8.4f}{ph:>8.4f}", file=out)
        else:
            print(f"\n  Perfect classification!", file=out)
        summary.append(("Likelihood (fine)", c, t))

    return summary


# ═══════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="fastlckin evaluation tool — compare output against ground truth",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single result
  python3 tests/evaluate.py tests/data/ground_truth.tsv tests/data/kinship_validation_output.tsv

  # Multiple results for comparison
  python3 tests/evaluate.py gt.tsv plink_out.tsv vcf_out.tsv

  # With custom labels
  python3 tests/evaluate.py --label PLINK --label VCF gt.tsv plink.tsv vcf.tsv
        """,
    )
    parser.add_argument("files", nargs="+",
                        help="ground_truth.tsv followed by one or more fastlckin output files")
    parser.add_argument("--label", action="append", default=[],
                        help="Label for each result file (optional, repeated per result)")
    parser.add_argument("--no-ibd", action="store_true",
                        help="Skip IBD accuracy report")
    parser.add_argument("--no-classify", action="store_true",
                        help="Skip classification report")
    args = parser.parse_args()

    if len(args.files) < 2:
        parser.error("Need at least ground_truth.tsv and one result file")

    gt_path = args.files[0]
    result_paths = args.files[1:]

    if not os.path.isfile(gt_path):
        print(f"{RED}ERROR: {gt_path} not found{RESET}")
        sys.exit(1)

    # Labels: use provided, or derive from filename
    labels = args.label if args.label else [os.path.splitext(os.path.basename(p))[0] for p in result_paths]
    while len(labels) < len(result_paths):
        labels.append(f"Result{len(labels)+1}")

    # Load data
    gt = load_ground_truth(gt_path)
    all_results = []
    for path in result_paths:
        if not os.path.isfile(path):
            print(f"{YELLOW}WARNING: {path} not found, skipping{RESET}")
            continue
        all_results.append((labels[len(all_results)], load_results(path)))

    print(f"Ground truth: {len(gt)} pairs")
    for lbl, res in all_results:
        print(f"  {lbl}: {len(res)} pairs from {os.path.basename(result_paths[all_results.index((lbl, res))])}")

    # ── Reports ───────────────────────────────────────────────────────
    grand_summary = []

    for lbl, results in all_results:
        print(f"\n{'╔' + '═'*68 + '╗'}")
        print(f"{'║'}  {lbl}".ljust(69) + f"{'║'}")
        print(f"{'╚' + '═'*68 + '╝'}")

        # IBD accuracy
        if not args.no_ibd:
            n_fail, n_miss = report_ibd_accuracy(gt, results, lbl)
            total = sum(1 for k in gt if k in results)
            n_pass = total - n_fail - n_miss
            grand_summary.append((f"{lbl} IBD", n_pass, total))

        # PI_HAT stats
        report_pihat_stats(gt, results, lbl)

        # Classification
        if not args.no_classify:
            cls_summary = report_classification(gt, results, lbl)
            for name, c, t in cls_summary:
                grand_summary.append((f"{lbl} {name}", c, t))

    # ── Grand Summary ─────────────────────────────────────────────────
    if len(all_results) > 1 or grand_summary:
        print(f"\n{'╔' + '═'*68 + '╗'}")
        print(f"{'║'}  GRAND SUMMARY".ljust(69) + f"{'║'}")
        print(f"{'╚' + '═'*68 + '╝'}")
        print(f"\n{'Method':<40}{'Accuracy':>14}{'Errors':>10}")
        print(f"{'-'*64}")
        for name, correct, total in grand_summary:
            pct = correct / total * 100 if total else 0
            print(f"{name:<40}{correct}/{total} ({pct:.1f}%){total - correct:>10}")


if __name__ == "__main__":
    main()
