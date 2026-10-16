#!/usr/bin/env python3
"""Run current submission modules against an authorized local input tree."""
import argparse
import json
import os
from pathlib import Path
import subprocess

BASE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', choices=['dm', 'matched'])
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--check-inputs', action='store_true')
    args = parser.parse_args()
    root = args.data_root.resolve()
    out = args.output.resolve()
    if not root.is_dir():
        parser.error('The data root must be an existing directory.')
    paths = [root / s for s in json.loads((BASE / 'required_inputs.json').read_text())[args.module]]
    repo = BASE.parent
    if (out == root or out == repo or repo in out.parents
            or any(p.parent == out or p.parent in out.parents for p in paths)
            or any(t in out.parts for t in ['inputs_private_20261007', 'private_materialized'])):
        parser.error('Output must be outside the code repository and source-input directories.')
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Output must be a new or empty directory.')
    missing = [str(s.relative_to(root)) for s in paths if not s.is_file()]
    result = {'module': args.module, 'checked_files': len(paths),
              'missing': missing, 'statistical_models_run': False}
    if missing:
        print(json.dumps(result))
        raise SystemExit(1)
    env = dict(os.environ, DISCO_RESULTS_ROOT=str(root), DISCO_RUN_OUTPUT=str(out))
    for file in sorted((BASE / args.module).glob('*.R')):
        subprocess.run([args.rscript, '--vanilla', '-e',
                        'parse(commandArgs(TRUE)[1]);cat("Parsed\\n")', str(file)],
                       check=True, stdout=subprocess.DEVNULL, env=env)
    for relative in ['results', 'tables', 'logs', 'qa', 'output/pdf', 'staging', 'input']:
        (out / relative).mkdir(parents=True, exist_ok=True, mode=0o700)
    out.chmod(0o700)
    if args.check_inputs:
        (out / 'results/input_check.json').write_text(json.dumps(result, indent=2))
        print(json.dumps(result))
        return
    names = ['fit_logDM.R', 'plot_Fig3c.R'] if args.module == 'dm' else ['fit_matched.R', 'plot_matched.R']
    for name in names:
        with (out / 'logs' / f'{name}.log').open('w') as log:
            subprocess.run([args.rscript, '--vanilla', str(BASE / args.module / name)],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(out)


if __name__ == '__main__':
    main()
