#!/usr/bin/env python3
"""Verify the reviewed public aggregate inputs in a checkout or source archive."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent


def verify(root=ROOT, require_tracked=False):
    root = Path(root).resolve()
    manifest = json.loads((root / 'PUBLIC_DATA_MANIFEST.json').read_text())
    tracked = None
    if require_tracked:
        tracked = set(subprocess.check_output(
            ['git', '-C', str(root), 'ls-files', '-z'], text=True).split('\0'))
    errors = []
    seen = set()
    for entry in manifest['files']:
        relative = entry['path']
        path = root / relative
        if relative in seen or Path(relative).is_absolute() or '..' in Path(relative).parts:
            errors.append(f'Invalid or duplicate manifest path: {relative}')
            continue
        seen.add(relative)
        if path.is_symlink() or not path.is_file() or root not in path.resolve().parents:
            errors.append(f'Missing or unsafe public input: {relative}')
            continue
        data = path.read_bytes()
        if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
            errors.append(f'Changed public input: {relative}')
        if tracked is not None and relative not in tracked:
            errors.append(f'Public input is not tracked by Git: {relative}')
    return errors, len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-tracked', action='store_true')
    args = parser.parse_args()
    errors, count = verify(require_tracked=args.require_tracked)
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {count} reviewed aggregate/contract files present with exact hashes')


if __name__ == '__main__':
    main()
