import gzip
import io
import math
import re
from functools import cached_property
from typing import Any, Literal

from pydantic import BaseModel
from pysam import AlignmentFile

# The regular expression pattern for "chrXX 1234567" in both VarScan2 Output and
# VCF files:
CONTIG_POSITION_PATTERN = re.compile(
    r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|[XY]|MT?)\t[0-9]+\b"
)
# More lenient pattern:
PATTERN_CHR_POSITION = re.compile(r"[^\t]+\t[0-9]+\b")
PATTERN_CHROM = re.compile(r"(?:chr)?([1-9]|1[0-9]|2[0-2]|[XY]|MT?)\W")

# Valid Phred+33 quality strings:
VALID_QUALITY_CHARS = [chr(33 + i) for i in range(42)]

# Define which chromosome coordinate is ahead for the following function:
CHROMOSOMES = [str(i) for i in range(1, 23)]
CHROMOSOMES.append("X")
CHROMOSOMES.append("Y")
CHROMOSOMES.append("M")


nan = float("nan")
inf = float("inf")

AA_3to1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}
AA_1to3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}


class VCFVariantRecord(BaseModel):
    chromosome: str | None = None
    position: int | None = None
    identifier: str | None = None
    refbase: str | None = None
    altbase: str | None = None
    qual: float | int | None = None
    filters: str | None = None
    info: str | None = None
    field: str | None = None
    samples: list[str] | None = None

    def get_info_items(self) -> list[str]:
        assert self.info
        return self.info.split(";")

    def get_info_value(self, variable: str) -> str | bool:
        if not self.info:
            raise ValueError("INFO field is empty.")
        key_item: Any = re.search(rf"\b{variable}=([^;\s]+)([;\W]|$)", self.info)
        # If key has a value attached to it, e.g., VAR=1,2,3, will return 1,2,3.
        if key_item:
            return key_item.groups()[0]
        # Perhaps it's simply a flag without "="
        else:
            key_item = self.info.split(";")
            return True if variable in key_item else False

    def get_sample_variable(self) -> list[str]:
        return self.field.split(":")

    def get_sample_item(
        self, idx: int = 0, out_type: Literal["dict", "list"] = "dict"
    ) -> dict[str, str] | tuple[list[str], list[str]]:
        """d to output a dictionary. l to output a tuple of lists"""
        assert self.samples
        if out_type.lower() == "dict":
            return dict(zip(self.get_sample_variable(), self.samples[idx].split(":")))
        elif out_type.lower() == "list":
            return (self.get_sample_variable(), self.samples[idx].split(":"))

    def get_sample_value(self, variable: str, idx: int = 0) -> dict[str, str]:
        assert self.field
        assert self.samples
        var2value = dict(zip(self.field.split(":"), self.samples[idx].split(":")))
        if variable not in var2value:
            return None
        return var2value[variable]

    @classmethod
    def from_vcf_line(cls, vcf_line: str) -> BaseModel:  # typing.Self in python>=3.11
        vcf_line = vcf_line.rstrip("\n")
        if not vcf_line:
            return cls()
        item = vcf_line.split("\t")
        (
            chromosome,
            pos,
            identifier,
            refbase,
            altbase,
            qual_str,
            filters,
            info,
            *has_samples,
        ) = item
        position = int(pos)
        if qual_str != ".":
            qual = float(qual_str)
        else:
            qual = None
        try:
            field, *samples = has_samples
        except ValueError:
            field = None
            samples = None
        return cls(
            chromosome=chromosome,
            position=position,
            identifier=identifier,
            refbase=refbase,
            altbase=altbase,
            qual=qual,
            filters=filters,
            info=info,
            field=field,
            samples=samples,
        )


class PysamHeader(BaseModel):
    bam_file: str

    @cached_property
    def bam_header(self) -> dict[str, Any]:
        bam = AlignmentFile(self.bam_file)
        return bam.header.to_dict()

    @cached_property
    def sample_name(self) -> tuple[str]:
        name_set = set()
        for header_i in self.bam_header["RG"]:
            name_set.add(header_i["SM"])
        name_tuple = tuple(name_set)
        return name_tuple  # type: ignore


def skip_vcf_header(opened_file):
    line_i = opened_file.readline().rstrip()
    while line_i.startswith("#"):
        line_i = opened_file.readline().rstrip()
    return line_i


def faiordict2contigorder(
    file_name: str, file_format: Literal["fai", "dict"]
) -> dict[str, int]:
    """Takes either a .fai or .dict file, and return a contig order dictionary,
    i.e., chrom_seq['chr1'] == 0"""

    if file_format not in ("fai", "dict"):
        raise TypeError("file_format has to be either fai or dict.")

    contig_sequence = []
    with open(file_name) as gfile:
        for line_i in gfile:
            if file_format == "fai":
                contig_match = re.match(r"([^\t]+)\t", line_i)
            elif file_format == "dict":
                if line_i.startswith("@SQ"):
                    contig_match = re.match(r"@SQ\tSN:([^\t]+)\tLN:", line_i)

            if contig_match:
                contig_i = contig_match.groups()[0].split(" ")[
                    0
                ]  # some .fai files have space after the contig.
                contig_sequence.append(contig_i)

    chrom_seq = {}
    for n, contig_i in enumerate(contig_sequence):
        chrom_seq[contig_i] = n

    return chrom_seq


def open_textfile(file_name: str) -> io.TextIOWrapper:
    # See if the input file is a .gz file:
    if str(file_name).lower().endswith(".gz"):
        return gzip.open(file_name, "rt")

    else:
        return open(file_name)


def open_bam_file(file_name: str) -> AlignmentFile | io.TextIOWrapper:
    try:
        return AlignmentFile(file_name, "rb")
    except ValueError:
        return open(file_name)


def ascii2phred33(x: str) -> int:
    """Put in an ASCII string, return a Phred+33 score."""
    if len(x) != 1:
        raise ValueError("Input is a single character.")
    return ord(x) - 33


def phred33toascii(x: int) -> str:
    """Put in a Phred33 score, return the character."""
    return chr(x + 33)


def p2phred(p, max_phred=inf):
    """Convert p-value to Phred-scale quality score."""
    if p < 0 or p > 1:
        raise ValueError("p-value must be between 0 and 1.")
    if p == 0:
        return max_phred
    if math.isnan(p):
        return nan
    if p == 1:
        return 0
    Q = -10 * math.log10(p)
    if Q > max_phred:
        Q = max_phred
    return Q


def phred2p(phred):
    """Convert Phred-scale quality score to p-value."""
    return 10 ** (-phred / 10)


def findall_index(mylist: list[Any], tolookfor: Any) -> list[int]:
    """Find all instances in a list that matches exactly thestring."""
    all_indices = [i for i, item in enumerate(mylist) if item == tolookfor]
    return all_indices


def findall_index_regex(mylist, pattern):
    """Find all instances in a list that matches a regex pattern."""
    all_indices = [i for i, item in enumerate(mylist) if re.search(pattern, item)]
    return all_indices


def count_repeating_bases(sequence):
    """For a string, count the number of characters that appears in a row.
    E.g., for string "ABBCCCDDDDAAAAAAA", the function returns 1, 2, 3, 4, 7,
    because there is 1 A, 2 B's, 3 C's, 4 D's, and then 7 A's.
    """
    counters = []
    previous_base = None
    for current_base in sequence:
        if current_base == previous_base:
            counters[-1] += 1
        else:
            counters.append(1)
        previous_base = current_base
    return counters


chrom_seq = {}
for n, contig_i in enumerate(CHROMOSOMES):
    chrom_seq[contig_i] = n


def whoisbehind(coord_0, coord_1, CHROMOSOMES):
    """
    coord_0 and coord_1 are two strings or two lists, specifying the chromosome,
    a (typically) tab, and then the location. Return the index where the
    coordinate is behind. Return 10 if they are the same position.
    """
    end_of_0 = end_of_1 = False
    if coord_0 == "" or coord_0 == ["", ""] or coord_0 == ("", "") or not coord_0:
        end_of_0 = True
    if coord_1 == "" or coord_1 == ["", ""] or coord_1 == ("", "") or not coord_1:
        end_of_1 = True
    if end_of_0 and end_of_1:
        return 10
    elif end_of_1:
        return 0
    elif end_of_0:
        return 1
    else:
        if isinstance(coord_0, str):
            chrom0, position0 = coord_0.split()
        elif isinstance(coord_0, list) or isinstance(coord_0, tuple):
            chrom0, position0 = coord_0[0], coord_0[1]
        if isinstance(coord_1, str):
            chrom1, position1 = coord_1.split()
        elif isinstance(coord_1, list) or isinstance(coord_1, tuple):
            chrom1, position1 = coord_1[0], coord_1[1]
        if isinstance(CHROMOSOMES, dict):
            chrom0_position = CHROMOSOMES[chrom0]
            chrom1_position = CHROMOSOMES[chrom1]
        elif isinstance(CHROMOSOMES, list) or isinstance(CHROMOSOMES, tuple):
            chrom0_position = CHROMOSOMES.index(chrom0)
            chrom1_position = CHROMOSOMES.index(chrom1)
        if chrom0_position < chrom1_position:
            return 0  # 1st coordinate is ahead
        elif chrom0_position > chrom1_position:
            return 1  # 1st coordinate is ahead

        # Must be in the same chromosome
        else:
            position0 = int(position0)
            position1 = int(position1)
            if position0 < position1:
                return 0
            elif position0 > position1:
                return 1
            # Same chromosome, same position, then same coordinate:
            elif position0 == position1:
                return 10


def vcf_header_modifier(infile_handle, addons=[], getlost=" "):
    """addons = A list of INFO, FORMAT, ID, or Filter lines you want to add.
    getlost = a regex expression for the ID of INFO/FORMAT/FILTER that you want
    to get rid of.
    """

    line_i = infile_handle.readline().rstrip()

    # First, write into the INFO and FORMAT what I want to add:
    vcfheader_info_format_filter = []
    vcfheader_misc = []
    for additions in addons:
        vcfheader_info_format_filter.append(additions)

    while line_i.startswith("##"):
        if re.match(r"##fileformat=", line_i):
            vcffileformat = line_i
        elif re.match(r"##(INFO|FORMAT|FILTER)", line_i):
            if not re.match(rf"##(INFO|FORMAT|FILTER)=<ID={getlost},", line_i):
                vcfheader_info_format_filter.append(line_i)
        elif re.match(r"##", line_i):
            vcfheader_misc.append(line_i)
        # Continue:
        line_i = infile_handle.readline().rstrip()

    # Print headers:
    vcfheader_info_format_filter.sort()
    vcfheader_misc.sort()
    return vcffileformat, vcfheader_info_format_filter, vcfheader_misc, line_i


def catchup(coordinate_i, line_j, filehandle_j, CHROMOSOMES):
    """
    Keep reading the j_th vcf file until it hits (or goes past) the i_th
    coordinate, at which time the function stops reading and you can do stuff.
    Returns (True, Vcf_line_j)  if the j_th vcf file contains an entry that
    matches the i_th coordinate. Returns (False, Vcf_line_j) if the j_th vcf
    file does not contain such an entry, and therefore the function has run past
    the i_th coordinate, by which time the programmer can decide to move into
    the next i_th coordiate.
    """
    coordinate_j = re.match(PATTERN_CHR_POSITION, line_j)
    if coordinate_j:
        coordinate_j = coordinate_j.group()
    else:
        coordinate_j = ""

    # Which coordinate is behind?
    is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

    # The file_j is already ahead, return the same line_j, but tag it "False"
    if is_behind == 0:
        reporter = (False, line_j)

    # The two coordinates are the same, return the same line_j, but tag it "True"
    elif is_behind == 10:
        reporter = (True, line_j)

    # If file_j is behind, then needs to catch up:
    elif is_behind == 1:
        # Keep at it until line_j is no longer behind:
        while is_behind == 1:
            # Catch up
            line_j = filehandle_j.readline().rstrip()
            next_coord = re.match(PATTERN_CHR_POSITION, line_j)
            if next_coord:
                coordinate_j = next_coord.group()
            else:
                coordinate_j = ""
            is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

        # If file_j has caught up exactly to the position of coordinate_i:
        if is_behind == 10:
            reporter = (True, line_j)

        # If file_j has run past coordinate_i:
        elif is_behind == 0:
            reporter = (False, line_j)

    return reporter


def catchup_multilines(coordinate_i, line_j, filehandle_j, CHROMOSOMES):
    """
    Keep reading the j_th vcf file until it hits (or goes past) the i_th coordinate, then
        1) Create a list to store information for this coordinate in the j_th
           vcf file
        2) Keep reading the j_th vcf file and store all lines with the same
           coordinate, until the coordinate goes to the next coordiate at which
           time the function stops reading and you can do stuff with the list
           created above.
        3) Basically, it won't stop when vcf_j reaches the coordinate, but only
           stop when vcf_j has gone beyond the coordinate.

    Returns (True, [Vcf_lines], line_j) if the j_th vcf file contains an entry
    that matches the i_th coordinate.
    Returns (False, []        , line_j) if the j_th vcf file does not contain
    such an entry, and therefore the function has run past the i_th coordinate,
    by which time the programmer can decide to move into the next i_th
    coordiate.
    """

    coordinate_j = re.match(PATTERN_CHR_POSITION, line_j)

    if coordinate_j:
        coordinate_j = coordinate_j.group()
    else:
        coordinate_j = ""

    # Which coordinate is behind?
    is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

    # The file_j is already ahead, return the same line_j, but tag it "False"
    if is_behind == 0:
        reporter = (False, [], line_j)

    # The two coordinates are the same, return the same line_j, but tag it "True"
    elif is_behind == 10:
        # Create a list, initiated with the current line:
        lines_of_coordinate_i = [line_j]

        while is_behind == 10:
            line_j = filehandle_j.readline().rstrip()
            next_coord = re.match(PATTERN_CHR_POSITION, line_j)

            if next_coord:
                coordinate_k = next_coord.group()
                if whoisbehind(coordinate_j, coordinate_k, CHROMOSOMES) == 1:
                    raise Exception(
                        "{} does not seem to be properly sorted: {} then {}.".format(
                            filehandle_j.name, coordinate_j, coordinate_k
                        )
                    )
                coordinate_j = coordinate_k
            else:
                coordinate_j = ""

            is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

            # If the next line (still) has the same coordinate:
            if is_behind == 10:
                lines_of_coordinate_i.append(line_j)

        reporter = (True, lines_of_coordinate_i, line_j)

    # If file_j is behind, then needs to catch up:
    # This is an opportunity to check if the vcf_j file is properly sorted, by asserting current line cannot be "behind" a subsequent line
    elif is_behind == 1:
        # Keep at it until line_j is no longer behind:
        while is_behind == 1:
            # Catch up
            line_j = filehandle_j.readline().rstrip()
            next_coord = re.match(PATTERN_CHR_POSITION, line_j)

            if next_coord:
                coordinate_k = next_coord.group()
                if whoisbehind(coordinate_j, next_coord.group(), CHROMOSOMES) == 1:
                    raise Exception(
                        "{} does not seem to be properly sorted: {} then {}.".format(
                            filehandle_j.name, coordinate_j, coordinate_k
                        )
                    )
                coordinate_j = coordinate_k
            else:
                coordinate_j = ""

            is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

        # If file_j has caught up exactly to the position of coordinate_i:
        if is_behind == 10:
            # Create a list, initiated with the current line:
            lines_of_coordinate_i = [line_j]
            while is_behind == 10:
                line_j = filehandle_j.readline().rstrip()
                next_coord = re.match(PATTERN_CHR_POSITION, line_j)
                if next_coord:
                    coordinate_k = next_coord.group()
                    if whoisbehind(coordinate_j, coordinate_k, CHROMOSOMES) == 1:
                        raise Exception(
                            "{} does not seem to be properly sorted: {} then {}.".format(
                                filehandle_j.name, coordinate_j, coordinate_k
                            )
                        )
                    coordinate_j = coordinate_k
                else:
                    coordinate_j = ""

                is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)
                # If the next line (still) has the same coordinate:
                if is_behind == 10:
                    lines_of_coordinate_i.append(line_j)

            reporter = (True, lines_of_coordinate_i, line_j)

        elif is_behind == 0:
            reporter = (False, [], line_j)

    return reporter


def find_vcf_at_coordinate(my_coordinate, latest_vcf_line, vcf_file_handle, chrom_seq):
    """Best used in conjunction with catchup_multilines.
    Given the current coordinate, the latest vcf_line from a vcf file, and the
    vcf file handle, it will return all the VCF variants (as VCF objects) at the
    given coordinate as a dictionary, where the key is the ( (contig, position),
    ref_base_i, alt_base_i ). If there are two ALT bases in a given VCF line,
    the output dictionary will include two copies of this VCF object, with two
    different keys, each representing a different ALT base.
    """
    latest_vcf_run = catchup_multilines(
        my_coordinate, latest_vcf_line, vcf_file_handle, chrom_seq
    )
    latest_vcf_here = latest_vcf_run[1]

    vcf_variants = {}
    if latest_vcf_run[0]:
        for vcf_line_i in latest_vcf_here:
            vcf_i = VCFVariantRecord.from_vcf_line(vcf_line_i)

            # Some VCF files wrongly uses "/" to separate different ALT's
            altbases = re.split(r"[,/]", vcf_i.altbase)
            for alt_i in altbases:
                vcf_variants[
                    ((vcf_i.chromosome, vcf_i.position), vcf_i.refbase, alt_i)
                ] = vcf_i
            assert my_coordinate[1] == vcf_i.position

    latest_vcf_line = latest_vcf_run[-1]
    return latest_vcf_run[0], vcf_variants, latest_vcf_line


# Read the 2nd file (i.e., filehandle_j) one line down if it's behind the i_th coordinate:
def catchup_one_line_at_a_time(coordinate_i, line_j, filehandle_j, CHROMOSOMES):
    """
    A sister program of catch_up, the difference is that the j_th file will be
    read only once if the coordinate is behind i, so that it allows the
    programmer a chance to do something for coordinates that only occurs in j,
    whereas the catch_up function will keep reading until it gets to to gets
    past i, so the programmer has no chance to do anything for coordinates that
    occur only in j. Return (0, Vcf_line_j)  if the coordinate_j matches
    coordinate_i. Return (1, Vcf_line_j)  if the coordinate_j is ahead
    coordinate_i. Return (-1, Vcf_line_j) if the coordinate_j is behind of
    coordinate_i.
    """

    coordinate_j = re.match(PATTERN_CHR_POSITION, line_j)
    if coordinate_j:
        coordinate_j = coordinate_j.group()
    else:
        coordinate_j = ""

    # Which coordinate is behind?
    is_behind = whoisbehind(coordinate_i, coordinate_j, CHROMOSOMES)

    # The file_j is already ahead, return the same line_j, but tag it "False"
    if is_behind == 0:
        reporter = (1, line_j)

    # The two coordinates are the same, return the same line_j, but tag it "True"
    elif is_behind == 10:
        reporter = (0, line_j)

    # If file_j is behind, then needs to catch up:
    elif is_behind == 1:
        # Read one line into file_j:
        line_j_next = filehandle_j.readline().rstrip()
        re.match(PATTERN_CHR_POSITION, line_j_next)
        reporter = (-1, line_j_next)

    return reporter