import bz2
import gzip
import lzma
import re
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Pattern, Sequence, Tuple

FASTQ_EXTENSION = r"(\.(fq|fastq)(\.gz)?)"
FASTQ_PATTERN = re.compile(fr"(.*){FASTQ_EXTENSION}")
FASTQ_R1_PATTERN = re.compile(fr"(.*)((_(R?))|(.read))(1)(_(\d+))?{FASTQ_EXTENSION}")

GROUPED_FASTQ_COLOR = "\033[01;32m"
UNGROUPED_COLOR = "\033[01;31m"
NO_COLOR = "\033[00m"


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
        return "\n".join([self.read_id, self.seq, self.unused, self.qual])


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
        raise ValueError(f"Path did not match R1 FASTQ pattern: {file_path}")
    return FASTQ_R1_PATTERN.sub(r"\1", file_path.name)


# noinspection PyPep8Naming
def get_rN_fastq(file_path: Path, n: int) -> Path:
    """
    @param file_path:
    @param n:
    @return:
    """
    if not FASTQ_R1_PATTERN.match(file_path.name):
        raise ValueError(f"Path did not match R1 FASTQ pattern: {file_path}")
    new_filename = FASTQ_R1_PATTERN.sub(fr"\1\g<2>{n}\7\9", file_path.name)
    return file_path.with_name(new_filename)


def create_match_find_funcs(
    pattern: Pattern,
) -> Tuple[Callable[[Path], bool], Callable[[Path], bool], Callable[[Path], Iterable[Path]],]:
    """
    :param pattern: File name regex
    :return: A 3-tuple of functions:
     [0] accepts a Path, returns whether the path name matches the regex
     [1] accepts a Path, returns whether the path name matches the regex
         and the file exists on disk
     [2] accepts a Path representing a directory, recursively walks the
         directory tree and returns Paths of all files matching the regex
    """

    def matches_regex(path: Path) -> bool:
        return bool(pattern.match(path.name))

    def matches_regex_and_is_file(path: Path) -> bool:
        return matches_regex(path) and path.is_file()

    def find_all_matching_files(directory: Path) -> Iterable[Path]:
        yield from filter(matches_regex_and_is_file, directory.glob("**/*"))

    return matches_regex, matches_regex_and_is_file, find_all_matching_files


is_fastq_r1, is_fastq_r1_file, find_r1_fastq_files = create_match_find_funcs(FASTQ_R1_PATTERN)
is_fastq, is_fastq_file, find_fastq_files = create_match_find_funcs(FASTQ_PATTERN)


def find_grouped_fastq_files(directory: Path, n: int, verbose=True) -> Iterable[Sequence[Path]]:
    """
    :param directories:
    :param n: number of FASTQ files to find; returns R1 through R{n}
    :param verbose:
    :return: Iterable of Sequence[Path]s, with n FASTQ Paths in each inner sequence
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        fastq_files = [r1_fastq_file]
        fastq_files.extend(get_rN_fastq(r1_fastq_file, i) for i in range(2, n + 1))

        if all(fq.is_file() for fq in fastq_files):
            if verbose:
                print(GROUPED_FASTQ_COLOR + f"Found group of {n} FASTQ files:" + NO_COLOR)
                for fq in fastq_files:
                    print(f"\t{fq}")
            yield fastq_files
        else:
            if verbose:
                print(UNGROUPED_COLOR + "Found ungrouped FASTQ file(s):" + NO_COLOR)
                for fq in fastq_files:
                    if fq.is_file():
                        print(f"\t{r1_fastq_file}")


def collect_fastq_files_by_directory(directory: Path) -> Dict[Path, List[Path]]:
    """
    Walk `directory`, finding all FASTQ files. Group these by the directory
    the files are in, so we can create the same directory structure for the
    output of FastQC.

    :param directory: Path to directory containing FASTQ files
    :return: Mapping of *relative* directory names to lists of absolute
      FASTQ file paths. The relative directory names are used to stage
      output directories matching the same structure as the input directory
      tree; the FASTQ paths are passed directly to FastQC
    """
    files_by_directory = defaultdict(list)
    for fastq_file in find_fastq_files(directory):
        relative_fastq = fastq_file.relative_to(directory)
        files_by_directory[relative_fastq.parent].append(fastq_file)
    return files_by_directory
