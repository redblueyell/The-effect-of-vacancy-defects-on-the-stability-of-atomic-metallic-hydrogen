#!/usr/bin/env python3
"""Batch-generate vac-pure *-number.txt files for timesteps 0,1000,... (200 files).

For each timestep t in 0..199000 step 1000:
 - if `vac-pure.<t>.dump` exists, compute pairwise distances and write `vac-pure.<t>-number.txt`.
 - otherwise record missing.
After processing, call `aggregate_vac_counts.py` to produce summary and plot.
"""
from pathlib import Path
import sys
import importlib.util
import subprocess
from pathlib import Path as _Path

# dynamic import of compute_vac_distances from same directory
_this_dir = _Path(__file__).parent
_cvd_path = _this_dir / 'compute_vac_distances.py'
spec = importlib.util.spec_from_file_location('compute_vac_distances', str(_cvd_path))
cvd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cvd)


def process_timestep(t, cwd: Path):
    dump_name = cwd / f'vac-pure.{t}.dump'
    out_name = cwd / f'vac-pure.{t}-number.txt'
    if not dump_name.exists():
        return False
    try:
        coords = cvd.parse_coords(dump_name)
        pairs = cvd.compute_pairwise(coords)
        isolated = cvd.atoms_all_greater_than(coords, 3.0)

        lines = []
        for a, b, d in pairs:
            lines.append(f'{a}-{b} {d:.2f}')
        lines.append('')
        if isolated:
            lines.append('Atoms with all distances > 3:')
            lines.append(' '.join(map(str, isolated)))
        else:
            lines.append('Atoms with all distances > 3: none')
        lines.append(f'Total count: {len(isolated)}')

        out_name.write_text('\n'.join(lines) + '\n')
        return True
    except Exception as e:
        print(f'Failed processing {dump_name}: {e}')
        return False


def main():
    cwd = Path('.')
    missing = []
    processed = []
    for i in range(200):
        t = i * 1000
        ok = process_timestep(t, cwd)
        if ok:
            processed.append(t)
            print(f'Processed timestep {t}')
        else:
            missing.append(t)

    print(f'Done. Processed {len(processed)} files, missing {len(missing)} files.')
    if missing:
        print('Missing timesteps (first 20):', missing[:20])

    # run aggregation
    try:
        subprocess.run([sys.executable, 'aggregate_vac_counts.py'], check=True)
        print('Aggregation completed.')
    except Exception as e:
        print('Aggregation failed:', e)


if __name__ == '__main__':
    main()
