# fastlckin

**Fast Maximum Likelihood Kinship Estimation from Low-Coverage Sequencing Data**

> A high-performance C++17 tool for maximum likelihood kinship estimation
> from low-coverage sequencing data using genotype likelihoods.

[![Build fastlckin](https://github.com/ShujiaHuang/fastlckin/actions/workflows/build.yml/badge.svg)](https://github.com/ShujiaHuang/fastlckin/actions/workflows/build.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

***fastlckin*** is a cross-platform, high-performance C++ tool for estimating
pairwise **kinship coefficients** (IBD sharing: k0, k1, k2) from **low-coverage
sequencing data** using genotype likelihoods. It implements the
Anderson & Weir (2007) IBS|IBD conditional probability model with FST
correction and uses Nelder-Mead maximum likelihood optimization — all without
requiring hard genotype calls.

> **Note on terminology:** "Genotype likelihood" in fastlckin refers to the
> general concept P(Data|Genotype), regardless of the VCF FORMAT field used.
> By default, the **PL** (Phred-scaled) field is preferred because it is more
> commonly available (e.g., GATK, DeepVariant). The GL (log10-scaled) field
> is used as a fallback. Use `--pl-field` to specify a custom field name.

fastlckin is built on top of [htslib](https://github.com/samtools/htslib)
(vendored as a git submodule) and uses multi-threaded processing via an
internal thread pool.

```bash
fastlckin v0.2.0
  fastlckin: Fast Maximum Likelihood Kinship Estimation from Low-Coverage Sequencing Data

Usage: fastlckin <command> [options]

Commands:
  relatedness   Estimate pairwise kinship (IBD coefficients)
  freq          Compute allele frequencies from PLINK .bed/.bim files
```

### Key improvements over lcMLkin v2.1

| Feature | lcMLkin (Python) | fastlckin (C++) |
| --- | --- | --- |
| Language | Python 3 + scipy | C++17 |
| Input modes | VCF + PLINK only | **VCF-only**, VCF + PLINK, **PLINK-only** |
| LD pruning | External PLINK call per pair | Built-in C++ (hard or expected genotypes) |
| VCF loading | `readlines()` loads all at once | Streaming via htslib (constant memory) |
| Optimization | `scipy.optimize.fmin` | Native Nelder-Mead (no Python dependency) |
| Parallelism | `multiprocessing` | Native `ThreadPool` |
| Temp files | Creates per-pair temp files | No temp files |
| Relationship classification | Not included | Built-in KING-based classification (`--classify`) |
| Bug fixes | Multiple known bugs (see spec) | All 13 known bugs resolved |

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

```bash
Usage: fastlckin <command> [options]

Commands:
  relatedness   Estimate pairwise kinship (IBD coefficients)
  freq          Compute allele frequencies from PLINK .bed/.bim files
```

---

## `fastlckin relatedness` — Pairwise kinship estimation

The main analysis command. It estimates IBD sharing coefficients (k0, k1, k2)
for every pair of individuals using the GLkin maximum likelihood framework.
Three input modes are supported:

### Input modes

| Mode | Input | AF source | LD pruning | Best for |
| --- | --- | --- | --- | --- |
| **VCF-only** | `-v <VCF>` | EM from GL | Expected genotypes from GL | Low-coverage sequencing, no reference panel |
| **VCF + PLINK** | `-v <VCF> -p <PREFIX>` | `.bed` hard genotypes | `.bed` hard genotypes | Low-coverage sequencing + reference panel |
| **PLINK-only** | `-p <PREFIX>` | `.bed` hard genotypes | `.bed` hard genotypes | SNP array / high-coverage data |

At least one of `-v` or `-p` must be provided. The mode is auto-detected
from the combination of arguments.

> **Mode selection guide:**
> - Low-coverage sequencing data without a reference panel → **VCF-only**
> - Low-coverage sequencing data with a reference panel → **VCF + PLINK**
> - SNP array or high-coverage data → **PLINK-only**
>
> All three modes produce the same output format and use the same
> relationship classification thresholds (see below).

### Full parameter reference of `relatedness`

```bash
Usage: fastlckin relatedness [-v <VCF>] [-p <PLINK_PREFIX>] [options]

At least one of -v or -p must be provided:
  -v, --vcf FILE          Input VCF/BCF file (.vcf, .vcf.gz, .bcf)
  -p, --plink PREFIX      PLINK binary file prefix (.bed/.bim/.fam)

Optional:
  -F, --freq FILE         Pre-computed .frq file (default: auto-compute)
  -f, --fst FLOAT         Prior FST value (default: 0.0)
  -t, --threads INT       Number of threads (default: 1)
  -o, --output FILE       Output TSV path (default: auto-generated)
      --maf-min FLOAT     Min allele frequency filter (default: 0.05)
      --maf-max FLOAT     Max allele frequency filter (default: 0.95)
      --ld-window INT     LD pruning window size in SNPs (default: 50)
      --ld-step INT       LD pruning step size (default: 5)
      --ld-r2 FLOAT       LD pruning r2 threshold (default: 0.5)
      --no-ld-prune       Skip LD pruning entirely (use all unmasked SNPs)
      --gq-min INT        Min GQ quality threshold (default: 1)
      --pl-field STR      VCF FORMAT field for Phred-scaled GL (default: PL)
      --n-restarts INT    Nelder-Mead restarts (default: 3)
      --xtol FLOAT        Optimizer parameter convergence (default: 0.01)
      --ftol FLOAT        Optimizer function convergence (default: 0.01)
      --classify          Enable automatic relationship classification
      --verbose           Verbose logging
  -h, --help              Show this help message
```

### Output format

The output is a TSV file with the following columns:

```tsv
Ind1	Ind2	k0	k1	k2	PI_HAT	N_SNPs	Relationship
SampleA	SampleB	1.000000	0.000000	0.000000	0.000000	12345	Unrelated
SampleA	SampleC	0.260000	0.490000	0.250000	0.495000	11987	First-degree
```

| Column | Description |
| --- | --- |
| `Ind1`, `Ind2` | Sample IDs |
| `k0` | IBD=0 probability |
| `k1` | IBD=1 probability |
| `k2` | IBD=2 probability |
| `PI_HAT` | Kinship coefficient: π̂ = 0.5·k1 + k2 |
| `N_SNPs` | Number of SNPs used (after LD pruning, or all unmasked SNPs with `--no-ld-prune`) |
| `Relationship` | Classification label (only with `--classify`) |

### Relationship classification (`--classify`)

fastlckin classifies pairwise relationships using PI_HAT (π̂ = 0.5·k1 + k2)
and IBD coefficients. The thresholds follow the log-scale midpoint scheme from
KING (Manichaikul et al., 2010), but adapted to the PI_HAT convention
(π̂ = 2φ, where φ is the kinship coefficient used by KING).

> These thresholds are **identical across all three input modes** (VCF-only,
> VCF+PLINK, PLINK-only) because PI_HAT is a property of the estimated IBD
> coefficients, not the input data type. However, different modes may produce
> slightly different PI_HAT estimates for the same pair due to differences
> in frequency estimation and LD pruning.

#### fastlckin classification rules

| PI_HAT range | Condition | Classification |
| --- | --- | --- |
| > 0.708 | k2 > 0.8 | `Duplicate/MZ` |
| [0.354, 0.708] | k0 < 0.05 | `Parent-Offspring` |
| [0.354, 0.708] | k0 ≥ 0.05 | `Full-Sibling` |
| [0.177, 0.354) | — | `Second-degree` |
| [0.0884, 0.177) | — | `Third-degree` |
| < 0.0884 | — | `Unrelated` |

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

### Usage examples

> Input data should be biallelic SNPs and autosomes only.

#### Mode 1 — VCF-only (no PLINK files needed)

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

#### Additional examples

**Using a pre-computed frequency file and custom FST (Mode 2):**

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -p /path/to/plink_prefix \
    -F allele_freq.frq \
    -f 0.01 \
    -t 8 \
    -o results.kinship.tsv
```

**Stricter LD pruning and MAF filtering:**

```bash
fastlckin relatedness \
    -v cohort.vcf.gz \
    -p /path/to/plink_prefix \
    --maf-min 0.05 --maf-max 0.95 \
    --ld-window 50 --ld-step 5 --ld-r2 0.2 \
    -t 8 \
    -o results.kinship.tsv
```

**Skipping LD pruning (for extremely sparse data like ancient DNA or NIPT):**

```bash
fastlckin relatedness \
    -v ancient_dna.vcf.gz \
    --no-ld-prune \
    -t 8 \
    -o results.kinship.tsv
```

When SNPs are very sparse, LD pruning removes valuable data and the r²
estimator itself is unreliable. `--no-ld-prune` retains all unmasked SNPs.

---

## `fastlckin freq` — Compute allele frequencies

Independently computes allele frequencies from PLINK binary files and
outputs a `.frq` file suitable for use with `relatedness -F`.

### Full parameter reference of `freq`

```bash
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

fastlckin implements the **GLkin** maximum likelihood framework:

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
   can be skipped entirely with `--no-ld-prune` (useful for extremely sparse
   data such as ancient DNA or NIPT where LD estimation is unreliable).

5. **Nelder-Mead Optimization** — Multi-restart Nelder-Mead simplex
   method finds the MLE of (k1, k2) under the IBD constraint space
   (k0 + k1 + k2 = 1, 4·k0·k2 ≤ k1²).

6. **Multi-threaded Processing** — All sample pairs are processed in
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
  structure).
- **Thread count**: Use `-t` to match available CPU cores — pair
  estimation is embarrassingly parallel.
- **GQ filtering**: `--gq-min 1` filters out truly missing genotypes.
  Increase to 10–20 (or even above) for stricter quality control on low-coverage data.

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

## Reference

- <https://genome.sph.umich.edu/w/images/1/1f/666.2017.12_-_Kinship_Coefficients.pdf>
- The original lcMLkin method: **Gazal et al. (2014)** *Bioinformatics*
- Anderson & Weir (2007) IBS|IBD model: **Anderson & Weir (2007)** *Genetics*
