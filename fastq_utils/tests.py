from pathlib import Path
import tempfile
from typing import List, Tuple
import unittest

from . import (
    collect_fastq_files_by_directory,
    is_fastq,
    is_fastq_r1,
    find_grouped_fastq_files,
    get_rN_fastq,
    get_sample_id_from_r1,
)

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
    ('prefix_0000.read1.fastq.gz', 'prefix_0000.read4.fastq.gz', 'prefix_0000'),
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

all_fastq_data_success_base = [
    'test.fq',
    'test.fq.gz',
    'test.fastq',
    'test.fastq.gz',
]
all_fastq_data_success_paths = [base_path / t for t in all_fastq_data_success_base]

all_fastq_data_failure_base = [
    'not_a_fastq.txt',
    'not_a_fq.zip',
]
all_fastq_data_failure_paths = [base_path / t for t in all_fastq_data_failure_base]

def touch(path: Path):
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, 'a'):
        pass

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

class TestIsFastq(unittest.TestCase):
    def test_success(self):
        for path in all_fastq_data_success_paths:
            with self.subTest(path=path):
                self.assertTrue(is_fastq(path))

    def test_failure(self):
        for path in all_fastq_data_failure_paths:
            with self.subTest(path=path):
                self.assertFalse(is_fastq(path))

class TestFindGroupedFastq(unittest.TestCase):
    def test_single_fastq(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            f = p / 'something_R1.fastq'
            touch(f)

            fastq_files = list(find_grouped_fastq_files(p, 1, verbose=False))
            self.assertEqual(len(fastq_files), 1)
            self.assertEqual(len(fastq_files[0]), 1)
            self.assertEqual(fastq_files[0][0], f)

    def test_four_fastqs(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            paths = [p / f'something_R{i}.fastq' for i in range(1, 5)]
            for path in paths:
                touch(path)
            touch(p / 'lone_R1.fastq')

            fastq_files = list(find_grouped_fastq_files(p, 4, verbose=False))
            self.assertEqual(len(fastq_files), 1)
            self.assertEqual(len(fastq_files[0]), 4)
            self.assertEqual(fastq_files[0], paths)

class TestCollectFastqByDirectory(unittest.TestCase):
    def test_collect_fastqs(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)

            subdir_count = 2
            file_count = 3
            dirs = [p / str(i) for i in range(subdir_count)]
            relative_dirs = sorted(d.relative_to(p) for d in dirs)

            for d in dirs:
                for j in range(file_count):
                    touch(d / f'{j}.fastq')

            grouped = collect_fastq_files_by_directory(p)
            self.assertEqual(sorted(grouped), relative_dirs)

            for relative_dir, files in grouped.items():
                with self.subTest(subdir=relative_dir):
                    self.assertEqual(len(files), file_count)
