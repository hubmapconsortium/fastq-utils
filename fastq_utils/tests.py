from pathlib import Path
from typing import List, Tuple
import unittest

from . import is_fastq_r1, get_rN_fastq, get_sample_id_from_r1

base_path = Path('path/to')
# there aren't any "R4" FASTQ files, just demonstrate generality
n = 4
# R1 FASTQ filename, R4 FASTQ filename, experiment base name
test_data_success_base: List[Tuple[str, str, str]] = [
    ('B001A001_1.fastq', 'B001A001_4.fastq', 'B001A001'),
    ('B001A001_1.fastq.gz', 'B001A001_4.fastq.gz', 'B001A001'),
    ('B001A001_1.fq', 'B001A001_4.fq', 'B001A001'),
    ('B001A001_1.fq.gz', 'B001A001_4.fq.gz', 'B001A001'),
    ('B001A001_R1.fastq', 'B001A001_R4.fastq', 'B001A001'),
    ('B001A001_R1.fastq.gz', 'B001A001_R4.fastq.gz', 'B001A001'),
    ('B001A001_R1.fq', 'B001A001_R4.fq', 'B001A001'),
    ('B001A001_R1.fq.gz', 'B001A001_R4.fq.gz', 'B001A001'),
    ('H4L1-4_S64_L001_R1_001.fastq.gz', 'H4L1-4_S64_L001_R4_001.fastq.gz', 'H4L1-4_S64_L001'),
]

def convert_success_data(t: Tuple[str, str, str]) -> Tuple[Path, Path, str]:
    return base_path / t[0], base_path / t[1], t[2]

test_data_success_paths = [convert_success_data(t) for t in test_data_success_base]

# not R1 FASTQ filenames
test_data_failure_base = [
    'H4L1-4_S64_L001_R2_001.fastq.gz',
    'B001A001_2.fq.gz',
]
test_data_failure_paths = [base_path / t for t in test_data_failure_base]

class TestIsFastqR1(unittest.TestCase):
    def test_success(self):
        for r1_path, r4_path, sample_id in test_data_success_paths:
            with self.subTest(r1_path=r1_path):
                self.assertTrue(is_fastq_r1(r1_path))

    def test_failure(self):
        for path in test_data_failure_paths:
            with self.subTest(path=path):
                self.assertFalse(is_fastq_r1(path))

class TestGetSampleID(unittest.TestCase):
    def test_success(self):
        for r1_path, r4_path, sample_id in test_data_success_paths:
            with self.subTest(r1_path=r1_path, sample_id=sample_id):
                self.assertEqual(get_sample_id_from_r1(r1_path), sample_id)

    def test_failure(self):
        for path in test_data_failure_paths:
            self.assertRaises(ValueError, get_sample_id_from_r1, path)

class TestGetRnFastq(unittest.TestCase):
    def test_success(self):
        for r1_path, r4_path, sample_id in test_data_success_paths:
            with self.subTest(r1_path=r1_path, r4_fastq=r4_path):
                self.assertEqual(get_rN_fastq(r1_path, n), r4_path)

    def test_failure(self):
        for path in test_data_failure_paths:
            self.assertRaises(ValueError, get_rN_fastq, path, n)
