#!/usr/bin/env python3
"""Convert PLINK binary files to VCF with PL fields, simulating low-coverage sequencing.

Usage: python3 plink_to_vcf_pl.py <prefix> <output_vcf> [--depth INT] [--seed INT]

Simulates sequencing reads at given average depth, computes genotype likelihoods
from allele counts, and outputs PL (phred-scaled likelihoods) in VCF FORMAT field.
"""

import struct, sys, math, random, argparse

def binomial(n, p, rng):
    """Simple binomial sampler using inverse CDF."""
    if p <= 0: return 0
    if p >= 1: return n
    # For small n, just sum Bernoulli trials
    return sum(1 for _ in range(n) if rng.random() < p)

def read_bed_genotypes(bed_path, n_samples, n_snps):
    """Read PLINK .bed file and return genotype matrix [snp][sample]."""
    # PLINK .bed format: 3-byte header, then packed genotypes
    # Each SNP: ceil(n_samples/4) bytes, 2 bits per sample
    # Bits: 00=hom_ref(0), 01=missing, 10=het(1), 11=hom_alt(2)
    genotypes = []
    bytes_per_snp = (n_samples + 3) // 4

    with open(bed_path, 'rb') as f:
        header = f.read(3)  # magic number + mode
        if header[:2] != b'\x6c\x1b':
            raise ValueError("Not a valid PLINK .bed file")
        mode = header[2]  # 0x01 = SNP-major

        for _ in range(n_snps):
            data = f.read(bytes_per_snp)
            if len(data) < bytes_per_snp:
                break
            snp_gt = []
            for i in range(n_samples):
                byte_idx = i // 4
                bit_idx = (i % 4) * 2
                bits = (data[byte_idx] >> bit_idx) & 0x03
                if bits == 0:
                    snp_gt.append(0)   # hom ref
                elif bits == 2:
                    snp_gt.append(1)   # het
                elif bits == 3:
                    snp_gt.append(2)   # hom alt
                else:
                    snp_gt.append(-1)  # missing
            genotypes.append(snp_gt)
    return genotypes

def read_bim(bim_path):
    """Read .bim file: chrom, snp_id, cm, pos, a1, a2."""
    snps = []
    with open(bim_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            snps.append({
                'chrom': parts[0], 'id': parts[1], 'cm': parts[2],
                'pos': int(parts[3]), 'a1': parts[4], 'a2': parts[5]
            })
    return snps

def read_fam(fam_path):
    """Read .fam file: FID, IID, PID, MID, sex, phenotype."""
    samples = []
    with open(fam_path) as f:
        for line in f:
            parts = line.strip().split()
            samples.append({'fid': parts[0], 'iid': parts[1]})
    return samples

def simulate_pl(true_gt, depth, rng):
    """Simulate PL from true genotype and sequencing depth.

    Models allele counts as Binomial(depth, allele_freq) where allele_freq
    depends on true genotype. Adds realistic noise for low-coverage data.
    """
    # True allele frequency given genotype
    if true_gt == 0:
        true_af = 0.02  # small error rate for hom ref
    elif true_gt == 1:
        true_af = 0.5
    else:
        true_af = 0.98  # small error rate for hom alt

    # Simulate read counts
    n_ref = binomial(depth, 1.0 - true_af, rng)
    n_alt = depth - n_ref

    # Compute genotype likelihoods (log-scale)
    # P(reads | GT) for GT in {0,1,2}
    # GT=0: all reads should be ref → P = (1-e)^n_ref * e^n_alt
    # GT=1: half ref, half alt → P = 0.5^depth
    # GT=2: all reads should be alt → P = e^n_ref * (1-e)^n_alt
    e = 0.01  # base error rate

    def log_likelihood(gt, n_ref, n_alt, depth, e):
        if gt == 0:
            # Expect all ref alleles
            return n_alt * math.log(e + 1e-10) + n_ref * math.log(1 - e + 1e-10)
        elif gt == 1:
            # Expect 50:50
            return depth * math.log(0.5)
        else:
            # Expect all alt alleles
            return n_ref * math.log(e + 1e-10) + n_alt * math.log(1 - e + 1e-10)

    ll = [log_likelihood(gt, n_ref, n_alt, depth, e) for gt in range(3)]

    # Convert to PL (phred-scaled, normalized so min=0)
    # PL = -10 * log10(P(GT|reads) / P(best_GT|reads))
    max_ll = max(ll)
    pls = []
    for l in ll:
        pl = int(round(-10 * (l - max_ll) / math.log(10)))
        pls.append(min(pl, 999))  # cap at 999

    return pls

def main():
    parser = argparse.ArgumentParser(description='Convert PLINK to VCF+PL')
    parser.add_argument('prefix', help='PLINK file prefix')
    parser.add_argument('output', help='Output VCF path')
    parser.add_argument('--depth', type=int, default=8, help='Average sequencing depth')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    args = parser.parse_args()

    rng = random.Random(args.seed)

    print(f"Reading PLINK files: {args.prefix}")
    samples = read_fam(f"{args.prefix}.fam")
    snps = read_bim(f"{args.prefix}.bim")
    n_samples = len(samples)
    n_snps = len(snps)
    print(f"  {n_samples} samples, {n_snps} SNPs")

    print("Reading genotypes from .bed...")
    genotypes = read_bed_genotypes(f"{args.prefix}.bed", n_samples, n_snps)
    print(f"  Read {len(genotypes)} SNPs")

    print(f"Simulating sequencing data (depth={args.depth})...")
    print(f"Writing VCF to: {args.output}")

    with open(args.output, 'w') as vcf:
        # VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=plink_to_vcf_pl\n")
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">\n')
        vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        vcf.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')

        # Header line
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples:
            vcf.write(f"\t{s['iid']}")
        vcf.write("\n")

        # Variants
        for snp_idx in range(n_snps):
            snp = snps[snp_idx]
            gts = genotypes[snp_idx]
            ref = snp['a2']  # PLINK a2 is usually the reference allele
            alt = snp['a1']

            vcf.write(f"{snp['chrom']}\t{snp['pos']}\t{snp['id']}\t{ref}\t{alt}\t.\tPASS\t.\tGT:PL:DP:GQ")

            for sample_idx in range(n_samples):
                true_gt = gts[sample_idx]
                if true_gt < 0:
                    vcf.write("\t./.:.:0:.")
                    continue

                # Simulate depth with some variation (Poisson-like)
                depth = max(1, int(rng.gauss(args.depth, args.depth * 0.3)))
                pls = simulate_pl(true_gt, depth, rng)

                # GT field
                if true_gt == 0:
                    gt_str = "0/0"
                elif true_gt == 1:
                    gt_str = "0/1"
                else:
                    gt_str = "1/1"

                # GQ = min non-zero PL
                gq = min(p for p in pls if p > 0) if any(p > 0 for p in pls) else 0
                pl_str = ",".join(str(p) for p in pls)

                vcf.write(f"\t{gt_str}:{pl_str}:{depth}:{gq}")

            vcf.write("\n")

            if (snp_idx + 1) % 1000 == 0:
                print(f"  {snp_idx + 1}/{n_snps} SNPs written")

    print(f"Done. {n_snps} SNPs x {n_samples} samples → {args.output}")

if __name__ == '__main__':
    main()
