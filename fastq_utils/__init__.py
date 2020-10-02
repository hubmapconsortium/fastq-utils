import bz2
from dataclasses import dataclass
from enum import Enum
import gzip
import lzma
from os import PathLike
from pathlib import Path
import re
from typing import Iterable, Sequence

FASTQ_EXTENSION = r'(\.(fq|fastq)(\.gz)?)'
FASTQ_PATTERN = re.compile(fr'(.*){FASTQ_EXTENSION}')
FASTQ_R1_PATTERN = re.compile(fr'(.*)_(R?)(1)(_(\d+))?{FASTQ_EXTENSION}')

GROUPED_FASTQ_COLOR = '\033[01;32m'
UNGROUPED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

# TODO: consolidate test data, make unit tests instead of doctests

class FileType(Enum):
    def __new__(cls, filetype, open_function):
        obj = object.__new__(cls)
        obj._value_ = filetype
        obj.open_function = open_function
        return obj

    GZ = ("gz", gzip.open)
    BZ2 = ("bz2", bz2.open)
    XZ = ("xz", lzma.open)
    TEXT = ("txt", open)

def get_file_type_by_extension(file_path: Path) -> FileType:
    suffix = file_path.suffix.lstrip(".")
    try:
        return FileType(suffix)
    except ValueError:
        # No special suffix, assume text
        return FileType.TEXT

def smart_open(file_path: PathLike, mode="rt", *args, **kwargs):
    file_type = get_file_type_by_extension(Path(file_path))
    return file_type.open_function(file_path, mode, *args, **kwargs)

@dataclass
class Read:
    read_id: str
    seq: str
    unused: str
    qual: str

    def serialize(self):
        return '\n'.join([self.read_id, self.seq, self.unused, self.qual])

revcomp_table = str.maketrans("ACTG", "TGAC")

def revcomp(seq: str) -> str:
    return seq.translate(revcomp_table)[::-1]

def fastq_reader(fastq_file: Path) -> Iterable[Read]:
    with smart_open(fastq_file) as f:
        while True:
            lines = [f.readline().strip() for _ in range(4)]
            if not all(lines):
                return
            yield Read(*lines)

def get_sample_id_from_r1(file_path: Path) -> str:
    """
    Only supports R1 FASTQ files.

    @param file_path:
    @return:
    """
    if not FASTQ_R1_PATTERN.match(file_path.name):
        raise ValueError(f'Path did not match R1 FASTQ pattern: {file_path}')
    return FASTQ_R1_PATTERN.sub(r'\1', file_path.name)

# noinspection PyPep8Naming
def get_rN_fastq(file_path: Path, n: int) -> Path:
    """
    @param file_path:
    @param n:
    @return:
    """
    if not FASTQ_R1_PATTERN.match(file_path.name):
        raise ValueError(f'Path did not match R1 FASTQ pattern: {file_path}')
    new_filename = FASTQ_R1_PATTERN.sub(fr'\1_\g<2>{n}\4\6', file_path.name)
    return file_path.with_name(new_filename)

def is_fastq_r1(path: Path) -> bool:
    """
    This is a separate (pure) function so it can have some unit tests
    :param path:
    :return: whether `fastq_file` is a R1 FASTQ file
    """
    return bool(FASTQ_R1_PATTERN.match(path.name))

def is_fastq_r1_file(path: Path) -> bool:
    """
    This is a separate (pure) function so it can have some unit tests
    :param path:
    :return: whether `fastq_file` is a R1 FASTQ file
    """
    return is_fastq_r1(path) and path.is_file()

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    yield from filter(is_fastq_r1_file, directory.glob('**/*'))

def is_fastq(path: Path):
    """
    :param path:
    :return: whether `path` is a FASTQ file
    """
    return bool(FASTQ_PATTERN.match(path.name))

def is_fastq_file(path: Path):
    return is_fastq(path) and path.is_file()

def find_all_fastq_files(directories: Iterable[Path]) -> Iterable[Path]:
    for directory in directories:
        yield from filter(is_fastq_r1_file, directory.glob('**/*'))

def find_grouped_fastq_files(directories: Iterable[Path], n: int, verbose=True) -> Iterable[Sequence[Path]]:
    """
    :param directories:
    :param n: number of FASTQ files to find; returns R1 through R{n}
    :param verbose:
    :return: Iterable of Sequence[Path]s, with n FASTQ Paths in each inner sequence
    """
    for directory in directories:
        for r1_fastq_file in find_r1_fastq_files(directory):
            fastq_files = [r1_fastq_file]
            fastq_files.extend(get_rN_fastq(r1_fastq_file, i) for i in range(2, n + 1))

            if all(fq.is_file() for fq in fastq_files):
                if verbose:
                    print(GROUPED_FASTQ_COLOR + f'Found group of {n} FASTQ files:' + NO_COLOR)
                    for fq in fastq_files:
                        print(f'\t{fq}')
                yield fastq_files
            else:
                if verbose:
                    print(UNGROUPED_COLOR + 'Found ungrouped FASTQ file(s):' + NO_COLOR)
                    for fq in fastq_files:
                        if fq.is_file():
                            print(f'\t{r1_fastq_file}')
