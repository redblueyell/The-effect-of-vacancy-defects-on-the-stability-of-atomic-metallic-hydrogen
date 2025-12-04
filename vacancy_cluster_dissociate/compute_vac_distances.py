#!/usr/bin/env python3
"""Compute pairwise distances for a vac-pure coordinate file.

Usage:
    python compute_vac_distances.py vac-pure.1000.dump

Output:
    Creates `<input>_distances.txt` with lines in format: `i-j distance` (distance with 2 decimals)
    Also prints and appends a list of atom indices whose distances to ALL other atoms are > 3 and the total count.
"""
import sys
import re
from pathlib import Path
from math import sqrt


def parse_coords(path: Path):
    """Parse coordinates from a file.

Supports simple formats like:
  - lines starting with index and three floats: `1 x y z`
  - LAMMPS dump ATOMS section (ITEM: ATOMS ...), will try to find x y z columns
"""
    text = path.read_text().splitlines()
    coords = []  # (index(int), x,y,z)
    box = None  # ( (xlo,xhi), (ylo,yhi), (zlo,zhi) ) or None

    # try simple index x y z lines first
    simple_re = re.compile(r'\s*(\d+)\s+([\-0-9\.eE]+)\s+([\-0-9\.eE]+)\s+([\-0-9\.eE]+)')
    for line in text:
        m = simple_re.match(line)
        if m:
            idx = int(m.group(1))
            x = float(m.group(2))
            y = float(m.group(3))
            z = float(m.group(4))
            coords.append((idx, x, y, z))

    if coords:
        # sort by provided index to get "按顺序编号"
        coords.sort(key=lambda t: t[0])
        # reindex to 1..N in file order
        reindexed = [(i+1, *coords[i][1:]) for i in range(len(coords))]
        return reindexed, box
    # attempt to parse BOX BOUNDS if present
    for i, line in enumerate(text):
        if line.startswith('ITEM: BOX BOUNDS'):
            lb = []
            j = i + 1
            while j < len(text) and len(lb) < 3:
                parts = text[j].split()
                if len(parts) >= 2:
                    try:
                        lo = float(parts[0]); hi = float(parts[1])
                        lb.append((lo, hi))
                    except Exception:
                        pass
                j += 1
            if len(lb) == 3:
                box = (lb[0], lb[1], lb[2])

    # fallback: try to parse LAMMPS dump ATOMS block (find header)
    atoms_idx = None
    for i, line in enumerate(text):
        if line.startswith('ITEM: ATOMS'):
            atoms_idx = i
            break
    if atoms_idx is None:
        raise ValueError(f'No coordinate lines found in {path}')

    header = text[atoms_idx].strip()
    cols = header.replace('ITEM: ATOMS', '').strip().split()
    # find x,y,z column indices
    try:
        xi = cols.index('x')
        yi = cols.index('y')
        zi = cols.index('z')
    except ValueError:
        raise ValueError('Could not find x,y,z columns in ATOMS header')

    # read following lines until next ITEM or end
    j = atoms_idx + 1
    while j < len(text) and not text[j].startswith('ITEM:'):
        toks = text[j].split()
        if len(toks) >= max(xi, yi, zi) + 1:
            # if id present at first column, use it else we assign sequential ids later
            try:
                idx = int(toks[0])
            except Exception:
                idx = None
            x = float(toks[xi])
            y = float(toks[yi])
            z = float(toks[zi])
            coords.append((idx if idx is not None else (j - atoms_idx), x, y, z))
        j += 1

    coords.sort(key=lambda t: t[0])
    reindexed = [(i+1, *coords[i][1:]) for i in range(len(coords))]
    return reindexed, box


def _min_image(dx, boxlen):
    if boxlen == 0:
        return dx
    return dx - round(dx / boxlen) * boxlen


def compute_pairwise(coords, box=None):
    n = len(coords)
    pairs = []
    if box is not None:
        lx = box[0][1] - box[0][0]
        ly = box[1][1] - box[1][0]
        lz = box[2][1] - box[2][0]
    else:
        lx = ly = lz = 0.0

    for i in range(n):
        idx_i, xi, yi, zi = coords[i]
        for j in range(i+1, n):
            idx_j, xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            if box is not None:
                dx = _min_image(dx, lx)
                dy = _min_image(dy, ly)
                dz = _min_image(dz, lz)
            d = sqrt(dx*dx + dy*dy + dz*dz)
            pairs.append((idx_i, idx_j, d))
    pairs.sort(key=lambda t: t[2])
    return pairs


def atoms_all_greater_than(coords, threshold, box=None):
    n = len(coords)
    res = []
    if box is not None:
        lx = box[0][1] - box[0][0]
        ly = box[1][1] - box[1][0]
        lz = box[2][1] - box[2][0]
    else:
        lx = ly = lz = 0.0

    for i in range(n):
        idx_i, xi, yi, zi = coords[i]
        ok = True
        for j in range(n):
            if i == j:
                continue
            _, xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            if box is not None:
                dx = _min_image(dx, lx)
                dy = _min_image(dy, ly)
                dz = _min_image(dz, lz)
            d = sqrt(dx*dx + dy*dy + dz*dz)
            if d <= threshold:
                ok = False
                break
        if ok:
            res.append(idx_i)
    return res


def main():
    if len(sys.argv) < 2:
        print('Usage: compute_vac_distances.py <vac-pure-file>')
        sys.exit(1)
    path = Path(sys.argv[1])
    if not path.exists():
        print(f'File not found: {path}')
        sys.exit(1)

    coords, box = parse_coords(path)
    pairs = compute_pairwise(coords, box)

    # allow optional output filename as second arg
    if len(sys.argv) >= 3:
        out_path = Path(sys.argv[2])
    else:
        out_path = path.with_name(path.stem + '_distances.txt')
    lines = []
    for a, b, d in pairs:
        # format: 编号-编号 距离 (两位小数)
        lines.append(f'{a}-{b} {d:.2f}')

    # atoms with all distances > 3
    isolated = atoms_all_greater_than(coords, 3.0, box)

    # append isolated info
    lines.append('')
    if isolated:
        lines.append('Atoms with all distances > 3:')
        lines.append(' '.join(map(str, isolated)))
    else:
        lines.append('Atoms with all distances > 3: none')
    lines.append(f'Total count: {len(isolated)}')

    out_path.write_text('\n'.join(lines) + '\n')
    print(f'Wrote {out_path}')
    print('First 40 lines:')
    for l in lines[:40]:
        print(l)


if __name__ == '__main__':
    main()
