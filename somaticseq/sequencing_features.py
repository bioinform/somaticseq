import pysam

import somaticseq.genomic_file_parsers.genomic_file_handlers as genome

nan = float("nan")


def get_homopolymer_lengths(
    ref_fa: pysam.FastaFile,
    my_coordinate: tuple[str, int, int],
    ref_base: str,
    first_alt: str,
) -> tuple[int, int]:
    """
    ref_fa is the opened reference fasta file handle
    my_coordiate is a list or tuple of 0-based (contig, position)
    """

    # Homopolymer eval (Make sure to modify for INDEL):
    # The min and max is to prevent the +/- 20 bases from exceeding the ends of
    # the reference sequence
    lseq = ref_fa.fetch(
        my_coordinate[0], max(0, my_coordinate[1] - 20), my_coordinate[1]
    )
    rseq = ref_fa.fetch(
        my_coordinate[0],
        my_coordinate[1] + 1,
        min(ref_fa.get_reference_length(my_coordinate[0]) + 1, my_coordinate[1] + 21),
    )
    # This is to get around some old version of pysam that reads the reference
    # sequence in bytes instead of strings
    lseq = (
        lseq.decode() if isinstance(lseq, bytes) else lseq  # type: ignore[attr-defined]
    )
    rseq = (
        rseq.decode() if isinstance(rseq, bytes) else rseq  # type: ignore[attr-defined]
    )
    seq41_ref = lseq + ref_base + rseq
    seq41_alt = lseq + first_alt + rseq
    ref_counts = genome.count_repeating_bases(seq41_ref)
    alt_counts = genome.count_repeating_bases(seq41_alt)
    homopolymer_length = max(max(ref_counts), max(alt_counts))

    # Homopolymer spanning the variant site:
    ref_c = 0
    alt_c = 0
    for i in rseq:
        if i == ref_base:
            ref_c += 1
        else:
            break

    for i in lseq[::-1]:
        if i == ref_base:
            ref_c += 1
        else:
            break

    for i in rseq:
        if i == first_alt:
            alt_c += 1
        else:
            break

    for i in lseq[::-1]:
        if i == first_alt:
            alt_c += 1
        else:
            break

    site_homopolymer_length = max(alt_c + 1, ref_c + 1)

    return homopolymer_length, site_homopolymer_length


def somatic_odds_ratio(
    n_ref: int, n_alt: int, t_ref: int, t_alt: int, max_value: float = 100
) -> float:
    # Odds Ratio just like VarDict's output
    sor_numerator = n_alt * t_ref
    sor_denominator = n_ref * t_alt
    if sor_numerator == 0 and sor_denominator == 0:
        sor = nan
    elif sor_denominator == 0:
        sor = max_value
    else:
        sor = sor_numerator / sor_denominator
        if sor >= max_value:
            sor = max_value
    return sor


def max_vocabularies(seq_length: int) -> int:
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    counts = 0
    k = 1
    while k <= seq_length:
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = int(
                counts + (seq_length - k + 1 + 1) * (seq_length - k + 1 - 1 + 1) / 2
            )
            break
        k += 1

    return counts


def linguistic_sequence_complexity(sequence: str) -> float:
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Assume 4 different nucleotides
    sequence = sequence.upper()
    if "N" in sequence:
        return float("nan")

    number_of_subseqs = 0
    seq_length = len(sequence)
    max_number_of_subseqs = max_vocabularies(seq_length)
    for i in range(1, seq_length + 1):
        set_of_seq_n = set()
        for n, nth_base in enumerate(sequence):
            if n + i <= len(sequence):
                sub_seq = sequence[n : n + i]
                set_of_seq_n.add(sub_seq)
        num_uniq_subseqs = len(set_of_seq_n)
        number_of_subseqs = number_of_subseqs + num_uniq_subseqs
    lc = number_of_subseqs / max_number_of_subseqs
    return lc


def max_sub_vocabularies(seq_length: int, max_subseq_length: int) -> int:
    # According to:
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # capping the length of sub_string as an input parameter
    assert max_subseq_length <= seq_length
    counts = 0
    k = 1
    while k <= max_subseq_length:
        if 4**k < (seq_length - k + 1):
            counts = counts + 4**k
        else:
            counts = int(
                counts
                + (2 * seq_length - k - max_subseq_length + 2)
                * (max_subseq_length - k + 1)
                / 2
            )
            break
        k += 1

    return counts


def ling_seq_complexity_with_max_vocab_length(
    sequence: str, max_substring_length: int = 20
) -> float:
    # Calculate linguistic sequence complexity according to
    # https://doi.org/10.1093/bioinformatics/18.5.679
    # Cut off substring at a fixed length
    sequence = sequence.upper()
    if "N" in sequence:
        return float("nan")

    number_of_subseqs = 0
    seq_length = len(sequence)
    max_number_of_subseqs = max_sub_vocabularies(seq_length, max_substring_length)
    set_of_seq_n: set[str] = set()
    for i in range(1, min(max_substring_length + 1, seq_length + 1)):
        set_of_seq_n.update(sequence[n : n + i] for n in range(len(sequence) - i + 1))
    number_of_subseqs = len(set_of_seq_n)
    lc = number_of_subseqs / max_number_of_subseqs
    return lc
