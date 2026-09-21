#!/usr/bin/env python3
"""
Extract 'e' values from NFS polynomial-selection output lines like:
    # norm 3.114939e-11 alpha -6.046667 e 3.713e-10 rroots 3
and plot a histogram of them with matplotlib.

Usage:
    python plot_e_histogram.py [input_file] [--bins N] [--log] [--out FILE.png]
"""
import re
import sys
import argparse
import matplotlib.pyplot as plt

# Matches the 'e <float>' token on '# norm ... e ... rroots ...' lines
E_PATTERN = re.compile(
    r'^#\s*norm\s+[-\d.eE+]+\s+alpha\s+[-\d.eE+]+\s+e\s+([-\d.eE+]+)\s+rroots\s+\d+'
)


def extract_e_values(path):
    values = []
    with open(path, "r") as f:
        for line in f:
            m = E_PATTERN.match(line.strip())
            if m:
                values.append(float(m.group(1)))
    return values


def main():
    ap = argparse.ArgumentParser(description="Plot histogram of NFS 'e' values.")
    ap.add_argument("input_file", nargs="?", default="nfs_dat.p",
                     help="Path to the poly-selection data file (default: nfs_dat.p)")
    ap.add_argument("--bins", type=int, default=50, help="Number of histogram bins")
    ap.add_argument("--log", action="store_true", help="Use log-scaled x-axis / bins")
    ap.add_argument("--out", default="e_histogram.png", help="Output image filename")
    args = ap.parse_args()

    values = extract_e_values(args.input_file)
    if not values:
        print("No 'e' values found — check the input file/format.", file=sys.stderr)
        sys.exit(1)

    print(f"Extracted {len(values)} e-values "
          f"(min={min(values):.3e}, max={max(values):.3e})")

    fig, ax = plt.subplots(figsize=(9, 6))

    if args.log:
        import numpy as np
        bins = np.logspace(np.log10(min(values)), np.log10(max(values)), args.bins)
        ax.set_xscale("log")
    else:
        bins = args.bins

    ax.hist(values, bins=bins, color="steelblue", edgecolor="black", alpha=0.85)
    ax.set_xlabel("e value")
    ax.set_ylabel("count")
    ax.set_title(f"Distribution of NFS polynomial 'e' scores (n={len(values)})")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved histogram to {args.out}")
    plt.show()


if __name__ == "__main__":
    main()
