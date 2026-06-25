# fastlckin

[![Build fastlckin](https://github.com/ShujiaHuang/fastlckin/actions/workflows/build.yml/badge.svg)](https://github.com/ShujiaHuang/fastlckin/actions/workflows/build.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

**Fast Maximum Likelihood Kinship Estimation from Low-Coverage Sequencing Data**

> A high-performance C++17 tool for maximum likelihood kinship estimation
> from low-coverage sequencing data using genotype likelihoods.

***fastlckin*** is a cross-platform, high-performance C++ tool for estimating
pairwise **kinship coefficients** (IBD sharing: k0, k1, k2) from **low-coverage
sequencing data** (0.1×–5×) using genotype likelihoods. It implements a
**dual-layer estimation framework**:

1. **Layer 1 — Continuous IBD estimation**: Anderson & Weir (2007) IBS|IBD
   conditional probability model with FST correction and Nelder-Mead MLE.
2. **Layer 2 — Exact model selection**: Per-SNP exact Mendelian likelihood
   ratios for 9 relationship types, aggregated by IBD equivalence classes.

All without requiring hard genotype calls.

> **Note on terminology:** "Genotype likelihood" in fastlckin refers to the
> general concept P(Data|Genotype), regardless of the VCF FORMAT field used.
> By default, the **PL** (Phred-scaled) field is preferred because it is more
> commonly available (e.g., GATK, DeepVariant). The GL (log10-scaled) field
> is used as a fallback. Use `--pl-field` to specify a custom field name.

fastlckin is built on top of [htslib](https://github.com/samtools/htslib)
(vendored as a git submodule) and uses multi-threaded processing via an
internal thread pool.

```
fastlckin v0.7.1
  fastlckin: Fast Maximum Likelihood Kinship Estimation from Low-Coverage Sequencing Data

Usage: fastlckin <command> [options]

Commands:
  relatedness   Estimate pairwise kinship (IBD coefficients)
  freq          Compute allele frequencies from PLINK .bed/.bim files
```

---

## Table of Contents

- [Installation](#installation)
- [Commands overview](#commands-overview)
- [fastlckin relatedness](#fastlckin-relatedness--pairwise-kinship-estimation)
- [fastlckin freq](#fastlckin-freq--compute-allele-frequencies)
- [Output format](#output-format)
- [Relationship classification](#relationship-classification)
- [Algorithm overview](#algorithm-overview)
- [Tips and best practices](#tips-and-best-practices)
- [Development](#development)
- [License](#license)
- [References](#references)

---

## Installation

### Option 1 — Download pre-built binary (Recommended, no compilation needed)

Pre-built static binaries are published on the
[GitHub Releases page](https://github.com/ShujiaHuang/fastlckin/releases).

| Platform              | Download                                                                                                          | Notes                                 |
| --------------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------- |
| Linux (x86_64)        | [fastlckin-linux-static](https://github.com/ShujiaHuang/fastlckin/releases/latest/download/fastlckin-linux-static) | Requires **glibc ≥ 2.35** (see below) |
| macOS (arm64 / Intel) | [fastlckin-macos-static](https://github.com/ShujiaHuang/fastlckin/releases/latest/download/fastlckin-macos-static) | Requires **macOS 12+**                |

#### System requirements for `fastlckin-linux-static`

The Linux binary is a **partial-static** build, produced on Ubuntu 22.04
(glibc 2.35) in CI. It bundles `libstdc++`, `libgcc`, `htslib`, `zlib`,
`bzip2`, and `xz` statically — only the system C library (**glibc**) is
linked dynamically. Because glibc symbol versions are **forward-compatible
only**, the binary requires the host `glibc` version to be **≥ 2.35**.

**Quick check on your machine:**

```bash
# If the printed glibc version is >= 2.35, fastlckin-linux-static will run.
ldd --version | head -1
```

```bash
# Linux
wget https://github.com/ShujiaHuang/fastlckin/releases/latest/download/fastlckin-linux-static
chmod +x fastlckin-linux-static
./fastlckin-linux-static --help
```

```bash
# macOS
curl -LO https://github.com/ShujiaHuang/fastlckin/releases/latest/download/fastlckin-macos-static
chmod +x fastlckin-macos-static
./fastlckin-macos-static --help
```

#### System requirements for `fastlckin-macos-static`

The macOS binary is a **best-effort static** build. Apple does not support
fully-static executables, so only Apple system frameworks (`libSystem`,
`libc++abi`, …) remain dynamic; everything else (htslib, zlib, bzip2, xz)
is statically linked.

- Tested on **macOS 13+** on Apple Silicon (arm64) — should also run on
  Intel macOS 12+.
- If you hit a "code signature invalid" error after `chmod +x`, run
  `xattr -d com.apple.quarantine ./fastlckin-macos-static` once to clear
  Gatekeeper's quarantine flag.

---

### Option 2 — Compile from source

*Requirements: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12,
autoconf/automake (for htslib), and the system libraries: zlib, bzip2,
xz-utils.*

#### Step 1 — Clone the repository (including the htslib submodule)

```bash
git clone --recursive https://github.com/ShujiaHuang/fastlckin.git
cd fastlckin
```

> If you forgot `--recursive`, run: `git submodule update --init --recursive`

#### Step 2 — Build with CMake (standard dynamic build)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable `bin/fastlckin` will be produced. Verify with:

```bash
./bin/fastlckin --help
./bin/fastlckin relatedness --help
```

#### Step 3 (Optional) — Run the unit tests

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/tests/generate_test_data   # generate synthetic test data
cd build && ctest --output-on-failure
```

#### Step 4 (Optional) — Build a static binary locally

**Linux** (portable static via Ubuntu/glibc — same approach used in CI):

```bash
sudo apt-get install -y build-essential cmake autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static --parallel
```

**macOS** (Homebrew):

```bash
brew install autoconf automake zlib bzip2 xz
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static --parallel
```

---

## Commands overview

```
Usage: fastlckin <command> [options]

Commands:
  relatedness   Estimate pairwise kinship (IBD coefficients)
  freq          Compute allele frequencies from PLINK .bed/.bim files

Use 'fastlckin <command> -h' for command-specific help.
```

---

## `fastlckin relatedness` — Pairwise kinship estimation

The main analysis command. It estimates IBD sharing coefficients (k0, k1, k2)
for every pair of individuals using a dual-layer maximum likelihood framework.

### Input modes

Three input modes are supported, auto-detected from the combination of `-v`
and `-p` arguments:

| Mode | Arguments | AF source | LD pruning | Best for |
| --- | --- | --- | --- | --- |
| **VCF-only** | `-v <VCF>` | EM from GL | Expected genotypes from GL | Low-coverage sequencing, no reference panel |
| **VCF + PLINK** | `-v <VCF> -p <PREFIX>` | `.bed` hard genotypes | `.bed` hard genotypes | Low-coverage sequencing + reference panel |
| **PLINK-only** | `-p <PREFIX>` | `.bed` hard genotypes | `.bed` hard genotypes | SNP array / high-coverage data |

At least one of `-v` or `-p` must be provided.

> **Mode selection guide:**
> - Low-coverage sequencing data without a reference panel → **VCF-only**
> - Low-coverage sequencing data with a reference panel → **VCF + PLINK**
> - SNP array or high-coverage data → **PLINK-only**

### Full parameter reference

```
Usage: fastlckin relatedness [-v <VCF>] [-p <PLINK_PREFIX>] [options]

At least one of -v or -p must be provided:
  -v, --vcf FILE          Input VCF/BCF file (.vcf, .vcf.gz, .bcf)
  -p, --plink PREFIX      PLINK binary file prefix (.bed/.bim/.fam)

Input files:
  -F, --freq FILE         Pre-computed .frq file (default: auto-compute)
      --pairs FILE        Estimate only specific pairs (TSV: ind1\tind2)
      --fst-file FILE     Per-SNP FST file (TSV: CHR\tPOS\tFST)

Population genetics:
  -f, --fst FLOAT         Prior FST value (default: 0.0)
      --maf-min FLOAT     Min allele frequency filter (default: 0.05)
      --maf-max FLOAT     Max allele frequency filter (default: 0.95)

LD pruning:
      --ld-window INT     LD pruning window size in SNPs (default: 50)
      --ld-step INT       LD pruning step size (default: 5)
      --ld-r2 FLOAT       LD pruning r² threshold (default: 0.5)
      --no-ld-prune       Skip LD pruning entirely (use all unmasked SNPs)

Quality control:
      --gq-min INT        Min GQ quality threshold (default: 1)
      --pl-field STR      VCF FORMAT field for Phred-scaled GL (default: PL)

Optimizer:
      --n-restarts INT    Nelder-Mead restarts (default: 5)
      --xtol FLOAT        Parameter convergence tolerance (default: 0.0001)
      --ftol FLOAT        Function convergence tolerance (default: 0.0001)

Relationship classification:
      --classify[=METHOD]  Classify relationship type (default method: threshold)
                           threshold: PI_HAT cutoffs (UN/PO/FS/2nd/3rd/DUP)
                           likelihood: per-SNP exact genotype likelihood ratio
      --classify-dup-threshold FLOAT     PI_HAT threshold for Duplicate/MZ (default: 0.708)
      --classify-first-threshold FLOAT   PI_HAT threshold for 1st degree (default: 0.354)
      --classify-second-threshold FLOAT  PI_HAT threshold for 2nd degree (default: 0.177)
      --classify-third-threshold FLOAT   PI_HAT threshold for 3rd degree (default: 0.0884)
      --output-likelihoods               Output per-relationship log-likelihood columns

KING-robust screening:
      --screen            Enable KING-robust screening (fast pre-filter)
      --screen-threshold FLOAT  PI_HAT threshold for screening (default: 0.0442)

Output:
  -o, --output FILE       Output TSV path (default: auto-generated)
  -t, --threads INT       Number of threads (default: 1)
      --verbose           Verbose logging
  -h, --help              Show this help message
```

### Usage examples

> Input data should be biallelic SNPs and autosomes only.

#### Quick start — VCF-only mode (simplest usage)

```bash
fastlckin relatedness \
    -v cohort_autosomes.vcf.gz \
    -t 8 \
    --classify \
    -o results.kinship.tsv
```

Allele frequencies are estimated from genotype likelihoods via the EM
algorithm. LD pruning uses posterior expected genotypes.

#### Mode 2 — VCF + PLINK (with reference panel)

```bash
fastlckin relatedness \
    -v cohort_autosomes.vcf.gz \
    -p /path/to/plink_prefix \
    -t 16 \
    --classify --verbose \
    -o results.kinship.tsv
```

Frequencies are computed from the PLINK `.bed` file (or from a pre-computed
`.frq` file via `-F`). LD pruning uses hard genotypes from `.bed`.

#### Mode 3 — PLINK-only (hard genotype mode)

```bash
fastlckin relatedness \
    -p /path/to/plink_prefix \
    -t 8 \
    --classify \
    -o results.kinship.tsv
```

Hard genotypes are treated as delta-function likelihoods in the same MLE
framework. Suitable for SNP array or high-coverage sequencing data.

#### Using likelihood-based classification (recommended for fine-grained relationships)

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -t 8 \
    --classify=likelihood \
    -o results.kinship.tsv
```

The `likelihood` method computes exact Mendelian likelihoods for 9
relationship types and selects the best one using IBD equivalence classes.
This provides finer-grained classification than `threshold` (see
[Relationship classification](#relationship-classification) below).

#### Using a pre-computed frequency file and custom FST

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -p /path/to/plink_prefix \
    -F allele_freq.frq \
    -f 0.01 \
    -t 8 \
    -o results.kinship.tsv
```

#### Estimating specific pairs only

```bash
# pairs.txt content (tab-separated):
# Ind1    Ind2
# SampleA SampleB
# SampleC SampleD

fastlckin relatedness \
    -v cohort.vcf.gz \
    --pairs pairs.txt \
    -t 4 \
    -o results.kinship.tsv
```

When analyzing a large cohort but only interested in specific pairs (e.g.,
known pedigrees), `--pairs` avoids the O(n²) all-pairs computation.

#### Using per-SNP FST values

```bash
# fst.txt content (tab-separated, with header):
# CHR  POS  FST
# chr1 10000 0.012
# chr1 20000 0.008

fastlckin relatedness \
    -v cohort.vcf.gz \
    --fst-file fst.txt \
    -t 8 \
    -o results.kinship.tsv
```

Per-SNP FST provides more granular population structure correction than a
single global `--fst` value.

#### KING-robust pre-screening (skip unrelated pairs)

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    --screen \
    --screen-threshold 0.0442 \
    -t 8 \
    --classify \
    -o results.kinship.tsv
```

`--screen` uses the KING-robust estimator as a fast pre-filter. Only pairs
with PI_HAT above the threshold (default: 0.0442, approximately 3rd-degree)
are analyzed with the full MLE. This dramatically reduces computation time
for large cohorts where most pairs are unrelated.

#### Stricter LD pruning and MAF filtering

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -p /path/to/plink_prefix \
    --maf-min 0.05 --maf-max 0.95 \
    --ld-window 50 --ld-step 5 --ld-r2 0.2 \
    -t 8 \
    -o results.kinship.tsv
```

#### Skipping LD pruning (for extremely sparse data)

```bash
fastlckin relatedness \
    -v ancient_dna.vcf.gz \
    --no-ld-prune \
    -t 8 \
    -o results.kinship.tsv
```

When SNPs are very sparse, LD pruning removes valuable data and the r²
estimator itself is unreliable. `--no-ld-prune` retains all unmasked SNPs.

#### Output per-relationship log-likelihoods

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -t 8 \
    --classify=likelihood \
    --output-likelihoods \
    -o results.kinship.tsv
```

This adds 9 additional columns (LL_UN, LL_DUP, LL_PO, LL_FS, LL_HS, LL_AV,
LL_GP, LL_FC, LL_GGP) with the log-likelihood for each relationship type.

---

## Output format

The output is a TSV file with comment headers and the following structure:

### Header comments (lines starting with `#`)

```tsv
# fastlckin relatedness v0.7.1
# Mode: VCF-only
# Command: fastlckin relatedness -v cohort.vcf.gz -t 8
# Date: 2026-06-24 10:30:00
# FST: 0  MAF_filter: [0.05, 0.95]  LD_prune: 50/5/0.5
```

### Column reference

| Column | Description |
| --- | --- |
| `Ind1`, `Ind2` | Sample IDs |
| `k0` | IBD=0 probability |
| `k1` | IBD=1 probability |
| `k2` | IBD=2 probability |
| `PI_HAT` | Kinship coefficient: π̂ = 0.5·k1 + k2 |
| `N_SNPs` | Number of SNPs used (after LD pruning, or all unmasked SNPs with `--no-ld-prune`) |
| `Relationship` | Classification label (with `--classify` or `--classify=likelihood`) |
| `MS_Relationship` | Model selection classification (with `--classify=likelihood`) |
| `LL_UN` ... `LL_GGP` | Per-relationship log-likelihoods (with `--output-likelihoods`) |
| `SE_k0`, `SE_k1`, `SE_k2`, `SE_PI_HAT` | Standard errors (from Fisher information matrix; NA if boundary) |
| `CI_k0_lo`, `CI_k0_hi` | 95% confidence interval for k0 |
| `CI_k1_lo`, `CI_k1_hi` | 95% confidence interval for k1 |
| `CI_k2_lo`, `CI_k2_hi` | 95% confidence interval for k2 |
| `CI_PI_lo`, `CI_PI_hi` | 95% confidence interval for PI_HAT (Delta method) |

### Example output

```tsv
Ind1	Ind2	k0	k1	k2	PI_HAT	N_SNPs	Relationship	SE_k0	SE_k1	SE_k2	SE_PI_HAT	CI_k0_lo	CI_k0_hi	CI_k1_lo	CI_k1_hi	CI_k2_lo	CI_k2_hi	CI_PI_lo	CI_PI_hi
SampleA	SampleB	1.000000	0.000000	0.000000	0.000000	12345	Unrelated	0.001200	NA	NA	0.000600	0.997648	1.002352	NA	NA	NA	NA	-0.001176	0.001176
SampleA	SampleC	0.000000	1.000000	0.000000	0.500000	11987	Parent-Offspring	NA	0.002100	NA	0.001050	NA	NA	0.995884	1.004116	NA	NA	0.497942	0.502058
SampleA	SampleD	0.250000	0.500000	0.250000	0.500000	11800	Full-Sibling	0.015000	0.020000	0.015000	0.012500	0.220601	0.279399	0.460801	0.539199	0.220601	0.279399	0.475500	0.524500
```

> **Note**: SE/CI values of `NA` indicate that the Fisher information matrix
> was singular or non-positive-definite at the estimate (common at boundary
> values like k0=0 or k2=0).

---

## Relationship classification

fastlckin supports two classification methods via `--classify[=METHOD]`:

### Method 1: Threshold-based classification (`--classify` or `--classify=threshold`)

Classifies relationships using PI_HAT (π̂ = 0.5·k1 + k2) and IBD
coefficients. The thresholds follow the log-scale midpoint scheme from
KING (Manichaikul et al., 2010).

| PI_HAT range | Condition | Classification |
| --- | --- | --- |
| > 0.708 | k2 > 0.8 | `Duplicate/MZ` |
| [0.354, 0.708] | k0 < 0.05 | `Parent-Offspring` |
| [0.354, 0.708] | k0 ≥ 0.05 | `Full-Sibling` |
| [0.177, 0.354) | — | `Second-degree` |
| [0.0884, 0.177) | — | `Third-degree` |
| < 0.0884 | — | `Unrelated` |

Custom thresholds can be set with `--classify-dup-threshold`,
`--classify-first-threshold`, `--classify-second-threshold`, and
`--classify-third-threshold`.

#### Comparison with KING thresholds

KING uses the **kinship coefficient** φ, while fastlckin reports **PI_HAT**
(π̂ = 0.5·k1 + k2 = 2φ). The numerical boundary values differ by a factor
of 2 because of this scale difference, but they correspond to the same
relationship categories.

| Relationship | fastlckin<br>PI_HAT (= 2φ) | KING<br>φ (kinship) | Expected 2φ | Expected φ | Additional criterion |
| --- | --- | --- | --- | --- | --- |
| Duplicate / MZ | > 0.708 | > 0.354 | 1.0 | 0.5 | fastlckin: k2 > 0.8; KING: π₀ < 0.1 |
| 1st degree (PO / FS) | [0.354, 0.708] | [0.177, 0.354] | 0.5 | 0.25 | PO vs FS: see below |
| 2nd degree | [0.177, 0.354) | [0.0884, 0.177) | 0.25 | 0.125 | — |
| 3rd degree | [0.0884, 0.177) | [0.0442, 0.0884) | 0.125 | 0.0625 | — |
| Unrelated | < 0.0884 | < 0.0442 | 0 | 0 | — |

#### Parent-Offspring vs Full-Sibling distinction

| | fastlckin | KING |
| --- | --- | --- |
| Metric | k0 (probability of IBD = 0) | π₀ (proportion of zero-IBS SNPs) |
| Parent-Offspring | k0 < 0.05 | π₀ < 0.1 |
| Full-Sibling | k0 ≥ 0.05 | π₀ ∈ (0.1, 0.365) |
| Rationale | Parent-offspring always share exactly one allele IBD at every locus (k0 = 0); full siblings may share 0, 1, or 2 alleles IBD (k0 ≈ 0.25) | Same principle using observed IBS = 0 proportion |

### Method 2: Likelihood-based classification (`--classify=likelihood`)

Computes exact Mendelian likelihoods for **9 relationship types** and
selects the best classification using **IBD equivalence classes**.

**Supported relationship types:**

| Relationship | IBD (k0, k1, k2) | Expected PI_HAT |
| --- | --- | --- |
| Unrelated | (1.0, 0.0, 0.0) | 0 |
| Duplicate/MZ | (0.0, 0.0, 1.0) | 1.0 |
| Parent-Offspring | (0.0, 1.0, 0.0) | 0.5 |
| Full-Sibling | (0.25, 0.5, 0.25) | 0.5 |
| Half-Sibling | (0.5, 0.5, 0.0) | 0.25 |
| Avuncular | (0.5, 0.5, 0.0) | 0.25 |
| Grandparent | (0.5, 0.5, 0.0) | 0.25 |
| First-Cousin | (0.75, 0.25, 0.0) | 0.125 |
| Great-Grandparent | (0.75, 0.25, 0.0) | 0.125 |

**IBD equivalence classes:** Under the "random spouse" assumption,
relationships with the same IBD coefficients produce **identical** genotype
joint probability distributions:

| Equivalence class | IBD (k0, k1, k2) | Contains |
| --- | --- | --- |
| Unrelated | (1.0, 0.0, 0.0) | Unrelated |
| Duplicate/MZ | (0.0, 0.0, 1.0) | Duplicate/MZ |
| Parent-Offspring | (0.0, 1.0, 0.0) | Parent-Offspring |
| Full-Sibling | (0.25, 0.5, 0.25) | Full-Sibling |
| 2nd-Degree | (0.5, 0.5, 0.0) | Half-Sibling, Grandparent, Avuncular |
| 3rd-Degree | (0.75, 0.25, 0.0) | First-Cousin, Great-Grandparent |

The likelihood method outputs two columns:
- `Relationship`: Coarse threshold-based classification (same as `threshold` method)
- `MS_Relationship`: Fine-grained likelihood-based classification using IBD
  equivalence classes (e.g., "2nd-Degree", "3rd-Degree")

> **When to use which method?**
> - `threshold`: Fast, simple, works well for clear-cut relationships.
>   Best for QC and quick overview.
> - `likelihood`: More informative, especially for 2nd/3rd-degree
>   relationships. Provides the statistically optimal classification.
>   Recommended for publication-quality results.

---

## `fastlckin freq` — Compute allele frequencies

Independently computes allele frequencies from PLINK binary files and
outputs a `.frq` file suitable for use with `relatedness -F`.

### Full parameter reference

```
Usage: fastlckin freq -p <PLINK_PREFIX> -o <OUTPUT.frq> [options]

Required:
  -p, --plink PREFIX   PLINK binary file prefix (.bed/.bim/.fam)
  -o, --output FILE    Output .frq file path

Optional:
  -t, --threads INT    Number of threads (default: 1)
  -h, --help           Show this help message
```

### Usage example

```bash
fastlckin freq \
    -p /path/to/plink_prefix \
    -o allele_freq.frq \
    -t 4
```

---

## Algorithm overview

fastlckin implements a **dual-layer estimation framework**:

### Layer 1: Continuous IBD estimation (IBS|IBD MLE)

1. **Genotype Likelihood Extraction** — Reads Phred-scaled genotype
   likelihoods from VCF (default field: PL; fallback: GL). The FORMAT field
   name is configurable via `--pl-field`. Values are converted to linear-scale
   likelihoods P(Data|G). In PLINK-only mode, hard genotypes are treated as
   delta-function likelihoods.

2. **Allele Frequency Estimation** — In VCF+PLINK and PLINK-only modes,
   frequencies are computed by counting alleles from `.bed` hard genotypes.
   In VCF-only mode, an **EM algorithm** estimates frequencies directly
   from genotype likelihoods using a Hardy-Weinberg prior:
   `P(G=0|p)=(1-p)²`, `P(G=1|p)=2p(1-p)`, `P(G=2|p)=p²`.

3. **Anderson & Weir IBS|IBD Model** — For each SNP, computes the
   conditional probability P(genotype_combo | IBD_state) using the
   Anderson & Weir (2007) formula with FST correction via the Mij function:
   `M(p, FST, i) = (1 - FST)·p + i·FST`.

4. **Built-in LD Pruning** — Sliding-window LD pruning (equivalent to
   PLINK `--indep-pairwise`) is implemented natively. In VCF+PLINK and
   PLINK-only modes, r² is computed from hard genotypes. In VCF-only mode,
   r² is computed from **posterior expected genotypes** E[G|D,p]. LD pruning
   can be skipped entirely with `--no-ld-prune`.

5. **Nelder-Mead Optimization** — Multi-restart Nelder-Mead simplex
   method finds the MLE of (k1, k2) under the IBD constraint space
   (k0 + k1 + k2 = 1, 4·k0·k2 ≤ k1²).

6. **Standard Errors** — Fisher information matrix (numerical Hessian)
   provides standard errors and 95% confidence intervals for k0, k1, k2,
   and PI_HAT (via Delta method).

### Layer 2: Exact model selection (v0.7.0)

7. **Exact Mendelian Likelihoods** — For each of the 9 relationship types,
   computes the exact genotype joint probability P(G1, G2 | R) using
   pedigree-specific Mendelian transmission formulas:

   | Relationship | Formula structure |
   | --- | --- |
   | Unrelated | Q(G1) × Q(G2) |
   | Duplicate/MZ | Q(G1) × δ(G1=G2) |
   | Parent-Offspring | 0.5 × Q(G1) × T(G2\|G1) + 0.5 × Q(G2) × T(G1\|G2) |
   | Full-Sibling | Σ Q(Gf) × Q(Gm) × T₂(G1\|Gf,Gm) × T₂(G2\|Gf,Gm) |
   | Half-Sibling | Σ Q(S) × T(G1\|S) × T(G2\|S) |
   | Grandparent | Q(G1) × Σ T(P\|G1) × T(G2\|P) |
   | Avuncular | Σ Q(GP) × T(G1\|GP) × [Σ T(P\|GP) × T(G2\|P)] |
   | First-Cousin | Σ Q(GP1) × Q(GP2) × [Σ T₂(P1\|GP1,GP2) × T(G1\|P1)] × [...] |
   | Great-Grandparent | Q(G1) × Σ T(A\|G1) × Σ T(B\|A) × T(G2\|B) |

   Where T = single-parent transmission, T₂ = two-parent transmission.

8. **IBD Equivalence Class Aggregation** — Log-likelihoods are aggregated
   by IBD equivalence classes, and the class with the highest likelihood
   is selected as the final classification.

### Common pipeline

9. **Multi-threaded Processing** — All sample pairs are processed in
   parallel using a native thread pool.

---

## Tips and best practices

- **PL over GL**: VCF files with PL fields are preferred — PL (Phred-scaled)
  is unambiguous and more commonly available (GATK, DeepVariant, etc.).
  GL (log10-scaled) is used as fallback. If your VCF uses a non-standard
  FORMAT field name for Phred-scaled genotype likelihoods, use
  `--pl-field <NAME>` to specify it.

- **Mode selection**: Use VCF-only mode when you only have a VCF. Use
  VCF+PLINK when a reference panel is available (more accurate frequencies).
  Use PLINK-only for SNP array data or high-coverage sequencing.

- **Classification method**: Use `--classify=likelihood` for
  publication-quality results. It provides statistically optimal
  classification using exact Mendelian likelihoods. Use `--classify`
  (threshold) for quick QC and overview.

- **KING-robust screening**: For large cohorts (>1000 samples), use
  `--screen` to skip the vast majority of unrelated pairs. This can
  reduce computation time by 10-100× with negligible loss of information.

- **Specific pairs**: When only interested in known pedigrees, use `--pairs`
  to avoid O(n²) all-pairs computation.

- **LD pruning**: The default parameters (`--ld-window 50 --ld-step 5
  --ld-r2 0.5`) are conservative. For dense SNP arrays or high-coverage
  (~30x or >20x) WGS data, consider tightening (`--ld-r2 0.2`). For extremely
  low-coverage data (e.g., ancient DNA, NIPT) where SNPs are very sparse,
  use `--no-ld-prune` to skip LD pruning entirely — the r² estimator is
  unreliable with too few overlapping samples, and pruning only discards
  valuable data.

- **VCF-only LD note**: In VCF-only mode, expected-genotype r² is a
  conservative (attenuated) estimator of true LD — the bias increases with
  lower coverage. This means fewer SNPs are pruned compared to hard-genotype
  LD pruning at the same `--ld-r2` threshold. For high-coverage VCF data,
  consider tightening (`--ld-r2 0.2–0.3`). For extremely low-coverage data,
  use `--no-ld-prune` instead — the attenuated r² estimator makes LD pruning
  largely ineffective when coverage is very low.

- **MAF filtering**: SNPs with very low or very high allele frequencies
  carry little information. The default `[0.05, 0.95]` range is suitable
  for most analyses.

- **FST correction**: For structured populations, set `--fst` to an
  appropriate value (e.g., 0.01 for mild structure, 0.05 for strong
  structure). For more granular correction, use `--fst-file` with per-SNP
  FST values.

- **Thread count**: Use `-t` to match available CPU cores — pair
  estimation is embarrassingly parallel.

- **GQ filtering**: `--gq-min 1` filters out truly missing genotypes.
  Increase to 10–20 (or even above) for stricter quality control on low-coverage data.

- **Standard errors**: SE values of `NA` in the output indicate that the
  estimate is at a boundary of the parameter space (e.g., k0=0 for
  Parent-Offspring). This is expected and does not indicate an error.

---

## Development

fastlckin is under active development. To update to the latest version:

```bash
git pull
git submodule update --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

To rebuild from scratch (recommended after a major htslib update):

```bash
rm -rf build
(cd htslib && make distclean) || true
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Pull requests, bug reports and feature requests are welcome at
<https://github.com/ShujiaHuang/fastlckin>.

---

## License

fastlckin is licensed under the
[GNU Affero General Public License v3.0](LICENSE).

---

## References

- Anderson, A. D., & Weir, B. S. (2007). A model-based approach for the
  estimation of relatedness using high-density genetic data. *Genetics*,
  175(2), 855–867.
- Balding, D. J., & Nichols, R. A. (1994). A method for quantifying
  differentiation between populations at multi-allelic loci and its
  implications for investigating identity and paternity. *Genetica*, 94(1), 3–12.
- Manichaikul, P., et al. (2010). Robust relationship inference in
  genome-wide association studies. *Bioinformatics*, 26(22), 2867–2873.
- Lipatov, M., et al. (2015). Maximum
  likelihood estimation of biological relatedness from low coverage
  sequencing data. *Genome Research*, 25(4), 459–466. (original lcMLkin method)
