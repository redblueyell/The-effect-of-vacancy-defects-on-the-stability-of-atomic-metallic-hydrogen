#!/usr/bin/env python3
"""Aggregate Total counts from vac-pure.*-number.txt and plot.

Outputs:
 - `vac-pure_counts_summary.txt` : human-readable summary (timestep -> count)
 - `vac-pure_counts.csv` : CSV with columns `timestep,count`
 - `vac-pure_counts.png` : plot (timestep vs count)

Usage:
    python aggregate_vac_counts.py
"""
from pathlib import Path
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def find_files(dirpath: Path):
    return sorted(dirpath.glob('vac-pure.*-number.txt'), key=lambda p: int(re.search(r'vac-pure\.(\d+)-number', p.name).group(1)) if re.search(r'vac-pure\.(\d+)-number', p.name) else p.name)


def extract_count(path: Path):
    text = path.read_text().splitlines()
    for line in reversed(text[-10:]):
        if 'Total count' in line:
            m = re.search(r'Total count:\s*(\d+)', line)
            if m:
                return int(m.group(1))
    # fallback: search entire file
    for line in text:
        if 'Total count' in line:
            m = re.search(r'Total count:\s*(\d+)', line)
            if m:
                return int(m.group(1))
    return None


def timestep_from_name(name: str):
    m = re.search(r'vac-pure\.(\d+)-number', name)
    if m:
        return int(m.group(1))
    return None


def main():
    cwd = Path('.')
    files = find_files(cwd)
    if not files:
        print('No vac-pure.*-number.txt files found in current directory.')
        sys.exit(1)

    data = []
    for p in files:
        t = timestep_from_name(p.name)
        c = extract_count(p)
        if c is None:
            print(f'Warning: could not find Total count in {p.name}')
            continue
        data.append((t, c))

    data.sort(key=lambda x: x[0])
    summary_txt = cwd / 'vac-pure_counts_summary.txt'
    csv_path = cwd / 'vac-pure_counts.csv'
    with summary_txt.open('w') as f:
        f.write('# Timestep -> Total count\n')
        for t, c in data:
            f.write(f'{t} -> {c}\n')

    with csv_path.open('w') as f:
        f.write('timestep,count\n')
        for t, c in data:
            f.write(f'{t},{c}\n')

    # plot
    ts = [t for t, _ in data]
    cs = [c for _, c in data]
    plt.figure(figsize=(10,4))
    plt.plot(ts, cs, marker='o', linestyle='-')
    plt.xlabel('Timestep')
    plt.ylabel('Total count (atoms with all distances > 3)')
    plt.title('vac-pure: atoms isolated by >3 distance')
    plt.grid(True)
    plt.tight_layout()
    out_png = cwd / 'vac-pure_counts.png'
    plt.savefig(out_png, dpi=150)
    print(f'Wrote {summary_txt}, {csv_path}, {out_png}')


if __name__ == '__main__':
    main()
