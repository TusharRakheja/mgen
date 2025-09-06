import argparse
import subprocess
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    prog='snplot',
    description='Generate genotype frequency heatmaps for populations and SNPs',
    epilog='Copyright (c) 2025 Tushar Rakheja (The MIT License)'
)

args = None

def parse_args():
    global args
    parser.add_argument('-i', '--input-set', required=True)
    parser.add_argument('-p', '--pop-labels', nargs='+', required=True)
    parser.add_argument('-s', '--snp-ids', nargs='+', required=True, help='Up to four SNP IDs')
    parser.add_argument('--plink', default="plink")
    parser.add_argument('--turn-on-wsl-for-plink', action='store_true')
    parser.add_argument('--savefig', type=str, default=None)
    args, _ = parser.parse_known_args()


def _normalize_list(tokens):
    out = []
    for token in tokens:
        out.extend([t for t in token.split(',') if t])
    return out


def _write_keep(fam_path, pop_id, keep_path):
    kept = 0
    with open(fam_path) as fam, open(keep_path, 'w') as keep:
        for line in fam:
            fields = line.strip().split()
            if fields[0] == pop_id:
                keep.write(f"{fields[0]} {fields[1]}\n")
                kept += 1
    return kept


def _run_plink(plink_cmd, use_wsl, bfile, snp, keep_file, out_prefix):
    cmd = (['wsl'] if use_wsl else []) + [
        plink_cmd,
        '--bfile', bfile,
        '--snp', snp,
        '--keep', keep_file,
        '--geno-counts',
        '--out', out_prefix
    ]
    subprocess.run(cmd, check=True)


def _read_counts(out_prefix):
    gcount = out_prefix + ".gcount"
    if not os.path.exists(gcount):
        raise FileNotFoundError(gcount)
    with open(gcount, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    header = lines[0][1:].split() if lines[0].startswith("#") else lines[0].split()
    data = lines[1].split()
    df = pd.DataFrame([data], columns=header)
    ref, alt = df.loc[0, ["REF", "ALT"]]
    counts = df.loc[0, ["HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS", "MISSING_CT"]].astype(int).tolist()
    labels = [f"{ref}{ref}", f"{ref}{alt}", f"{alt}{alt}", "Missing"]
    return labels, counts


def main():
    parse_args()
    plink_prefix = args.input_set
    pops = _normalize_list(args.pop_labels)
    snps = _normalize_list(args.snp_ids)[:4]

    if len(snps) > 4:
        print("Error: Maximum of four SNPs supported.", file=sys.stderr)
        sys.exit(1)

    fam_path = plink_prefix + ".fam"
    keep_file = "keep.txt"
    out_prefix = "plink_out"

    data = {}  # snp -> (xlabels, rows (percent), raw_counts, pops)
    for snp in snps:
        rows = []
        raw = []
        used_pops = []
        labels = None
        for pop in pops:
            kept = _write_keep(fam_path, pop, keep_file)
            if kept == 0:
                continue
            _run_plink(args.plink, args.turn_on_wsl_for_plink, plink_prefix, snp, keep_file, out_prefix)
            labels_, counts = _read_counts(out_prefix)
            total = sum(counts)
            if total == 0:
                continue
            rows.append([c/total for c in counts])
            raw.append(counts)
            used_pops.append(pop)
            labels = labels or labels_
            for ext in (".gcount", ".log"):
                path = out_prefix + ext
                if os.path.exists(path):
                    os.remove(path)
        data[snp] = (labels, rows, raw, used_pops)

    if os.path.exists(keep_file):
        os.remove(keep_file)

    n = len(snps)
    nrows, ncols = 2, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 10))
    colormaps = ['Blues', 'Greens', 'Reds', 'Purples']

    for idx, snp in enumerate(snps):
        r, c = divmod(idx, 2)
        ax = axes[r][c]
        labels, percents, raw, used_pops = data[snp]
        mat = np.array(percents)
        im = ax.imshow(mat, aspect='auto', cmap=colormaps[idx], vmin=0, vmax=1)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_yticks(range(len(used_pops)))
        ax.set_yticklabels(used_pops)
        ax.set_title(f"{snp}")
        ax.set_xlabel("Genotype")
        ax.set_ylabel("Population")
        for r in range(mat.shape[0]):
            for c in range(mat.shape[1]):
                ax.text(c, r, f"{mat[r,c]:.2f}\n({raw[r][c]})", ha='center', va='center', fontsize=9, color='white' if mat[r,c] > 0.70 else 'black')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="% of individuals")

    plt.tight_layout()
    if args.savefig:
        plt.savefig(args.savefig, dpi=200)
    else:
        plt.show()

if __name__ == "__main__":
    main()
