#!/usr/bin/env python3
# mypy: ignore-errors

import argparse

import pysam

parser = argparse.ArgumentParser(
    description=(
        "Count some metrics from BAM files such as "
        "fragment size, duplication rates, fraction of soft-clipped "
        "and discordant reads."
    ),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-bam",
    "--bam-file-in",
    type=str,
    help="Input BAM file",
    required=True,
    default=None,
)
parser.add_argument(
    "-maxl",
    "--max-length",
    type=int,
    help="max frag length to consider",
    required=False,
    default=1000,
)


args = parser.parse_args()
bam_file = args.bam_file_in
max_length = args.max_length

with pysam.AlignmentFile(bam_file) as bam:
    reads = bam.fetch()

    clipped_and_discordant = clipped_only = discordant_only = concordant_reads = mq0 = (
        unmapped
    ) = duplicated_reads = total_reads = 0
    frag_lengths = {}
    duplicates_per_length = {}
    MQs = {}

    for read_i in reads:
        if read_i.is_proper_pair:
            concordant_reads += 1
            frag_length = abs(read_i.template_length)

            if frag_length in frag_lengths:
                frag_lengths[frag_length] += 1

                if read_i.is_duplicate:
                    try:
                        duplicates_per_length[frag_length] += 1
                    except KeyError:
                        duplicates_per_length[frag_length] = 1

            else:
                frag_lengths[frag_length] = 1

                if read_i.is_duplicate:
                    duplicates_per_length[frag_length] = 1
                else:
                    duplicates_per_length[frag_length] = 0

        if read_i.is_duplicate:
            duplicated_reads += 1

        mq = read_i.mapping_quality
        if mq in MQs:
            MQs[mq] += 1
        else:
            MQs[mq] = 1

        if not read_i.is_unmapped:
            if (not read_i.is_proper_pair) and ("S" in read_i.cigarstring):
                clipped_and_discordant += 1
            elif not read_i.is_proper_pair:
                discordant_only += 1
            elif "S" in read_i.cigarstring:
                clipped_only += 1

            if read_i.mapping_quality == 0:
                mq0 += 1

        else:
            unmapped += 1

        total_reads += 1

    # Find fragment length median:
    n_reads_processed = 0
    for frag_i in sorted(frag_lengths):
        n_reads_processed += frag_lengths[frag_i]

        if n_reads_processed >= concordant_reads / 2:
            median_frag_length = frag_i
            break

    # Calculate mean fragment length
    total_length = 0
    total_reads_processed = 0
    for frag_i in frag_lengths:
        if 0 < frag_i < max_length:
            total_length += frag_i * frag_lengths[frag_i]
            total_reads_processed += frag_lengths[frag_i]

    mean_length = total_length / total_reads_processed

    # Calculate standard deviation of fragment length
    sum_of_square_of_x_minus_mean = 0
    for frag_i in frag_lengths:
        if 0 < frag_i < max_length:
            square_of_x_minus_mean = (frag_i - mean_length) ** 2
            sum_of_square_of_x_minus_mean += (
                square_of_x_minus_mean * frag_lengths[frag_i]
            )

    frag_length_std_dev = (sum_of_square_of_x_minus_mean / total_reads_processed) ** (
        1 / 2
    )

    print(
        "Duplicated reads and rates: {}, {}".format(
            duplicated_reads, duplicated_reads / total_reads
        )
    )
    print(f"soft-clipped and discordant reads: {clipped_and_discordant}")
    print(f"soft-clipped and concordant reads: {clipped_only}")
    print(f"discordant and not-clipped reads: {discordant_only}")
    print(f"MQ0 reads: {mq0}")
    print(f"unmapped reads: {unmapped}")
    print(f"Total reads: {total_reads}")
    print(f"Mean fragment length: {mean_length}")
    print(f"fragment length standard deviation: {frag_length_std_dev}")
    print(f"median fragment length: {median_frag_length}")

    print("###---\nMQ: Number, Fraction")
    for mq_i in sorted(MQs):
        print(f"MQ={mq_i}: {MQs[mq_i]}, {MQs[mq_i] / total_reads}")

    print("###---\nFrag length distribution:")
    for frag_i in sorted(frag_lengths):
        print(
            "FragLength={}: {} / {} = {}".format(
                frag_i,
                duplicates_per_length[frag_i],
                frag_lengths[frag_i],
                duplicates_per_length[frag_i] / frag_lengths[frag_i],
            )
        )
