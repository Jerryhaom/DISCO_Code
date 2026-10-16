"""Missing inputs and unsafe output paths must fail before models or writes."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'submission_20261013/run.py'


class SubmissionRunnerTests(unittest.TestCase):
    def invoke(self, module, data, output):
        return subprocess.run([sys.executable, str(RUNNER), module, '--data-root', str(data),
                               '--output', str(output), '--rscript', 'unavailable-Rscript'],
                              capture_output=True, text=True)

    def test_missing_inputs_stop_before_output_or_R(self):
        with tempfile.TemporaryDirectory(prefix='disco paths with spaces ') as tmp:
            data = Path(tmp) / 'data'; data.mkdir()
            for module, count in [('dm', 8), ('matched', 6)]:
                out = Path(tmp) / module
                result = self.invoke(module, data, out)
                self.assertEqual(result.returncode, 1)
                self.assertEqual(len(json.loads(result.stdout)['missing']), count)
                self.assertFalse(out.exists())
                self.assertNotIn('FileNotFoundError', result.stderr)

    def test_earlier_runners_reject_repository_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=ROOT/'unsafe-test-output'
            commands=[
                [sys.executable,str(ROOT/'analysis_20261008/run.py'),'theory-identities','--output',str(output)],
                [sys.executable,str(ROOT/'approved/run.py'),'ipaq','--data-root',tmp,'--output',str(output)]
            ]
            for command in commands:
                result=subprocess.run(command,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertFalse(output.exists())
                self.assertTrue('repository' in result.stderr or 'upload candidate' in result.stderr)

    def test_output_rejections(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / 'data'; data.mkdir()
            existing = Path(tmp) / 'existing'; existing.mkdir(); (existing / 'keep').write_text('preserve')
            cases = [data, ROOT / 'unsafe-test-output', data / 'private_materialized/new', existing]
            for output in cases:
                with self.subTest(output=output):
                    result = self.invoke('dm', data, output)
                    self.assertEqual(result.returncode, 2)
            self.assertEqual((existing / 'keep').read_text(), 'preserve')
            self.assertFalse((ROOT / 'unsafe-test-output').exists())


if __name__ == '__main__':
    unittest.main()
