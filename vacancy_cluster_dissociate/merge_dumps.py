#!/usr/bin/env python3
"""Merge voids-relax.*.dump and vac-pure.*.dump by timestep.

For each timestep present in both sets:
- Read `voids-relax.<t>.dump` and `vac-pure.<t>.dump`.
- Ensure merged `ITEM: ATOMS` header contains `type` column.
- Set type=1 for atoms from voids-relax, type=2 for atoms from vac-pure.
- Write `merged.<t>.dump` with updated NUMBER OF ATOMS and merged atom lines.

Usage: run in the directory containing the dump files:
    python merge_dumps.py
"""
import re
from pathlib import Path
import sys


def parse_dump(path: Path):
    text = path.read_text().splitlines()
    # capture useful blocks: TIMESTEP, NUMBER OF ATOMS, BOX BOUNDS (header + lines), ATOMS header
    timestep = None
    num_atoms = None
    box_header = None
    box_bounds = []
    atoms_header = None
    atoms_idx = None
    for i, line in enumerate(text):
        if line.strip() == 'ITEM: TIMESTEP':
            if i + 1 < len(text):
                timestep = text[i + 1].strip()
        if line.strip() == 'ITEM: NUMBER OF ATOMS':
            if i + 1 < len(text):
                try:
                    num_atoms = int(text[i + 1].strip())
                except Exception:
                    num_atoms = None
        if line.startswith('ITEM: BOX BOUNDS'):
            box_header = line.strip()
            j = i + 1
            while j < len(text) and not text[j].startswith('ITEM:') and len(box_bounds) < 10:
                box_bounds.append(text[j])
                j += 1
        if line.startswith('ITEM: ATOMS') and atoms_idx is None:
            atoms_header = line.strip()
            atoms_idx = i

    if atoms_idx is None:
        raise ValueError(f'No "ITEM: ATOMS" found in {path}')

    header_lines = text[:atoms_idx + 1]
    atom_lines = text[atoms_idx + 1:]
    # limit atom_lines to num_atoms if present
    if num_atoms is not None:
        atom_lines = atom_lines[:num_atoms]

    cols = []
    if atoms_header:
        m = re.match(r'ITEM: ATOMS\s*(.*)', atoms_header)
        if m:
            rest = m.group(1).strip()
            if rest:
                cols = rest.split()

    return {
        'path': path,
        'timestep': timestep,
        'header_lines': header_lines,
        'atom_lines': atom_lines,
        'num_atoms': num_atoms if num_atoms is not None else len(atom_lines),
        'box_header': box_header,
        'box_bounds': box_bounds,
        'cols': cols,
    }


def build_merged_columns(cols_voids, cols_vac):
    # Ensure the merged columns contain at least x,y,z and type; produce a canonical ordering
    merged_set = []
    for c in cols_voids:
        if c not in merged_set:
            merged_set.append(c)
    for c in cols_vac:
        if c not in merged_set:
            merged_set.append(c)

    # Ensure x,y,z present
    for coord in ('x', 'y', 'z'):
        if coord not in merged_set:
            merged_set.append(coord)
    # Ensure type present
    if 'type' not in merged_set:
        merged_set.append('type')

    # Prefer order: id, type, x, y, z, then others
    others = [c for c in merged_set if c not in ('id', 'type', 'x', 'y', 'z')]
    ordered = []
    ordered.append('id')
    if 'type' in merged_set:
        ordered.append('type')
    ordered.extend([c for c in ('x', 'y', 'z') if c in merged_set])
    ordered.extend(others)
    return ordered


def line_to_map(line, cols):
    toks = line.strip().split()
    m = {}
    for i, c in enumerate(cols):
        if i < len(toks):
            m[c] = toks[i]
        else:
            m[c] = '0'
    # store any extra tokens as numbered extras
    if len(toks) > len(cols):
        for j in range(len(cols), len(toks)):
            m[f'_extra{j-len(cols)}'] = toks[j]
    return m


def map_to_line(m, cols):
    toks = [m.get(c, '0') for c in cols]
    # append extras if present
    extras = [v for k, v in sorted(m.items()) if k.startswith('_extra')]
    toks.extend(extras)
    return ' '.join(toks)


def merge_for_timestep(t, voids_path: Path, vac_path: Path, out_path: Path):
    v = parse_dump(voids_path)
    p = parse_dump(vac_path)

    merged_cols = build_merged_columns(v['cols'], p['cols'])

    merged_atoms = []

    # voids atoms -> type=1
    for line in v['atom_lines']:
        m = line_to_map(line, v['cols'])
        m['type'] = '1'
        merged_atoms.append(map_to_line(m, merged_cols))

    # vac atoms -> type=2
    for line in p['atom_lines']:
        m = line_to_map(line, p['cols'])
        m['type'] = '2'
        merged_atoms.append(map_to_line(m, merged_cols))
    # Assign sequential unique ids for merged atoms (1..N)
    atom_maps = []
    # rebuild as maps to ensure consistent columns
    for line in v['atom_lines']:
        m = line_to_map(line, v['cols'])
        m['type'] = '1'
        atom_maps.append(m)
    for line in p['atom_lines']:
        m = line_to_map(line, p['cols'])
        m['type'] = '2'
        atom_maps.append(m)

    # assign ids
    for i, m in enumerate(atom_maps, start=1):
        m['id'] = str(i)

    merged_atom_lines = [map_to_line(m, merged_cols) for m in atom_maps]

    # build output header in canonical LAMMPS dump order
    out_lines = []
    # TIMESTEP
    timestep_val = v.get('timestep') or p.get('timestep') or str(t)
    out_lines.append('ITEM: TIMESTEP')
    out_lines.append(str(timestep_val))
    # NUMBER OF ATOMS
    out_lines.append('ITEM: NUMBER OF ATOMS')
    out_lines.append(str(len(merged_atom_lines)))
    # BOX BOUNDS: prefer header with flags if present
    if v.get('box_header'):
        out_lines.append(v['box_header'])
        out_lines.extend(v['box_bounds'])
    elif p.get('box_header'):
        out_lines.append(p['box_header'])
        out_lines.extend(p['box_bounds'])
    else:
        # fallback
        out_lines.append('ITEM: BOX BOUNDS')
        if v['box_bounds']:
            out_lines.extend(v['box_bounds'])
        elif p['box_bounds']:
            out_lines.extend(p['box_bounds'])

    # ATOMS header
    out_lines.append('ITEM: ATOMS ' + ' '.join(merged_cols))
    out_lines.extend(merged_atom_lines)

    out_path.write_text('\n'.join(out_lines) + '\n')


def main():
    cwd = Path('.')
    voids = {}
    vacs = {}
    for p in cwd.iterdir():
        if p.is_file():
            m = re.match(r'voids-relax\.(\d+)\.dump$', p.name)
            if m:
                voids[m.group(1)] = p
            m2 = re.match(r'vac-pure\.(\d+)\.dump$', p.name)
            if m2:
                vacs[m2.group(1)] = p

    common = sorted(set(voids.keys()) & set(vacs.keys()), key=lambda x: int(x))
    if not common:
        print('No matching timesteps found between voids-relax.*.dump and vac-pure.*.dump')
        return

    print(f'Found {len(common)} matching timesteps. Merging...')
    for t in common:
        out_name = f'merged.{t}.dump'
        print(f' Merging timestep {t} -> {out_name}')
        try:
            merge_for_timestep(t, voids[t], vacs[t], Path(out_name))
        except Exception as e:
            print(f'  Failed to merge {t}: {e}')

    print('Done.')


if __name__ == '__main__':
    main()
