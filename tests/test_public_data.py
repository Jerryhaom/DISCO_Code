"""Regression check for omitted/modified aggregate files in source archives."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('public_data', ROOT / 'verify_public_data.py')
public_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(public_data)


class PublicDataTests(unittest.TestCase):
    def test_release_has_all_aggregates(self):
        errors, count = public_data.verify()
        self.assertGreater(count, 6)
        self.assertEqual(errors, [])

    def test_missing_and_changed_input_rejected(self):
        entry = json.loads((ROOT / 'PUBLIC_DATA_MANIFEST.json').read_text())['files'][0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'PUBLIC_DATA_MANIFEST.json').write_text(json.dumps({'files': [entry]}))
            self.assertIn('Missing', public_data.verify(root)[0][0])
            path = root / entry['path']
            path.parent.mkdir(parents=True)
            shutil.copyfile(ROOT / entry['path'], path)
            self.assertEqual(public_data.verify(root)[0], [])
            path.write_bytes(path.read_bytes() + b'\n')
            self.assertIn('Changed', public_data.verify(root)[0][0])


if __name__ == '__main__':
    unittest.main()
