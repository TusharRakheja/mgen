import argparse
import subprocess
import pandas as pd
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    prog='snplot',
    description='Generate genotype frequency plots for populations',
    epilog='Copyright (c) 2025 Tushar Rakheja (The MIT License)'
)

args = None

def parse_args():
    global args
    parser.add_argument(
        '-i', '--input-set', dest='input_set', type=str, required=True,
        help='Prefix of your input set (e.g. just "prefix" if prefix.bed, prefix.bim, prefix.fam)'
    )
    parser.add_argument(
        '-p', '--pop-label', dest='pop_label', type=str, required=True,
        help='Population or family label'
    )
    parser.add_argument(
        '-s', '--snp-id', dest='snp_id', type=str, required=True,
        help='The SNP you want to plot the genotype frequencies for (e.g. rs12345)'
    )
    parser.add_argument('--plink', dest='plink', default="plink", type=str, help='Command used to invoke plink on your machine (e.g. "plink", "plink2" etc)')
    parser.add_argument('--turn-on-wsl-for-plink', default=False, action='store_true', help="(If running on Windows) Use Plink via WSL")
    args, _ = parser.parse_known_args()

def main():
    parse_args()

    plink_prefix = args.input_set
    snp_id = args.snp_id
    pop_id = args.pop_label

    keep_file = "keep.txt"
    out_prefix = "plink_out"

    # Write keep file from .fam
    fam_path = plink_prefix + ".fam"
    with open(fam_path, 'r') as fam, open(keep_file, 'w') as keep:
        for line in fam:
            line = line.strip()
            if not line:
                continue

            fields = line.split()
            if fields[0] == pop_id:
                keep.write(f"{fields[0]} {fields[1]}\n")

    # Run PLINK to get genotype counts
    subprocess.run([
        *(['wsl'] if args.turn_on_wsl_for_plink else []),
        args.plink,
        "--bfile", plink_prefix,
        "--snp", snp_id,
        "--keep", keep_file,
        "--geno-counts",
        "--out", out_prefix
    ], check=True)

    log_file = out_prefix + ".log"

    # Read counts
    counts_file = out_prefix + ".gcount"
    df = pd.read_csv(counts_file, delim_whitespace=True)

    # Extract allele names
    ref, alt = df.loc[0, ["REF", "ALT"]]

    # Build genotype labels
    genotypes = [
        f"{ref}{ref}",   # HOM_REF_CT
        f"{ref}{alt}",   # HET_REF_ALT_CTS
        f"{alt}{alt}",   # TWO_ALT_GENO_CTS
        "Missing"        # MISSING_CT
    ]
    counts = df.loc[0, ["HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS", "MISSING_CT"]].values

    # Plot
    plt.bar(genotypes, counts)
    plt.title(f"Genotype counts for {snp_id} in {pop_id}")
    plt.ylabel("Count")
    plt.xlabel("Genotype")
    plt.show()

    if os.path.exists(keep_file):
        os.remove(keep_file)

    if os.path.exists(counts_file):
        os.remove(counts_file)

    if os.path.exists(log_file):
        os.remove(log_file)

if __name__ == "__main__":
    main()
