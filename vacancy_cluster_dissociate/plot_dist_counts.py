#!/usr/bin/env python3
"""
plot_dist_counts.py

Reads `corrected_s1/dist_counts_per_1000.txt` (format: header then `step count` lines)
and plots `count` vs `step`. Saves `dist_counts_per_1000.png` into the same folder.

Usage:
  python plot_dist_counts.py --dir ./corrected_s1

Optional: --show to open the plot window.
"""
import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_summary(path):
    steps = []
    counts = []
    with open(path, 'r', encoding='utf-8') as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.lower().startswith('step'):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                step = int(parts[0])
                cnt = int(parts[1])
            except ValueError:
                continue
            steps.append(step)
            counts.append(cnt)
    return steps, counts


def plot(steps, counts, outpath, show=False):
    plt.figure(figsize=(10,4))
    plt.plot(steps, counts, marker='o', linestyle='-', color='tab:blue')
    plt.xlabel('Step')
    plt.ylabel('Number of distances > 3.0')
    plt.title('vac-pure: number of pairwise distances > 3.0 per 1000 steps')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    if show:
        # attempt to open in a window (won't work in headless CI)
        plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', default='corrected_s1', help='Directory containing dist_counts_per_1000.txt')
    parser.add_argument('--show', action='store_true', help='Show plot window')
    args = parser.parse_args()

    summary_path = os.path.join(args.dir, 'dist_counts_per_1000.txt')
    if not os.path.exists(summary_path):
        print(f"Summary file not found: {summary_path}")
        sys.exit(2)

    steps, counts = read_summary(summary_path)
    if not steps:
        print('No data parsed from summary file.')
        sys.exit(0)

    outpath = os.path.join(args.dir, 'dist_counts_per_1000.png')
    plot(steps, counts, outpath, show=args.show)
    print(f'Plot saved to: {outpath}')


if __name__ == '__main__':
    main()
