import math
from typing import Any, cast

import pysam

from somaticseq.bam_features import BamFeatures, collect_bam_features_batch


def assert_bam_features_equal(actual: BamFeatures, expected: BamFeatures) -> None:
    for field_name, expected_value in expected.model_dump().items():
        actual_value = getattr(actual, field_name)
        if isinstance(expected_value, float) and math.isnan(expected_value):
            assert math.isnan(actual_value)
        else:
            assert actual_value == expected_value


class RecordingAlignmentFile:
    def __init__(self, alignment_file: pysam.AlignmentFile) -> None:
        self.alignment_file = alignment_file
        self.fetch_calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, stop: int):
        self.fetch_calls.append((contig, start, stop))
        return self.alignment_file.fetch(contig, start, stop)


class PrefixedReadAlignmentFile:
    def __init__(self, prefix_reads: list[object], alignment_file: pysam.AlignmentFile) -> None:
        self.prefix_reads = prefix_reads
        self.alignment_file = alignment_file

    def fetch(self, contig: str, start: int, stop: int):
        for read in self.prefix_reads:
            yield read
        yield from self.alignment_file.fetch(contig, start, stop)


def test_collect_bam_features_batch_matches_per_variant_oracle_and_reduces_fetches(
    tiny_tumor_bam: str,
) -> None:
    with pysam.AlignmentFile(tiny_tumor_bam) as bam:
        recording_bam = RecordingAlignmentFile(bam)
        candidates = [
            (("1", 8450), "T", "G"),
            (("1", 8492), "G", "A"),
            (("1", 10387), "A", "C"),
            (("1", 5000), "A", "T"),
        ]

        batch_features = collect_bam_features_batch(cast(Any, recording_bam), candidates)

        assert len(recording_bam.fetch_calls) == 3
        for candidate in candidates:
            expected = BamFeatures.from_alignment_file(
                bam,
                candidate[0],
                candidate[1],
                candidate[2],
            )
            assert_bam_features_equal(
                batch_features[(candidate[0], candidate[1], candidate[2])],
                expected,
            )


def test_collect_bam_features_batch_skips_reads_without_reference_end(
    tiny_tumor_bam: str,
) -> None:
    class FakeRead:
        is_unmapped = False
        reference_start = 8440
        reference_end = None

    with pysam.AlignmentFile(tiny_tumor_bam) as bam:
        prefixed_bam = PrefixedReadAlignmentFile([FakeRead()], bam)
        candidate = (("1", 8450), "T", "G")

        batch_feature = collect_bam_features_batch(cast(Any, prefixed_bam), [candidate])[candidate]
        expected = BamFeatures.from_alignment_file(
            bam,
            candidate[0],
            candidate[1],
            candidate[2],
        )

        assert_bam_features_equal(batch_feature, expected)


def test_collect_bam_features_batch_handles_same_coordinate_alt_alleles_and_cluster_boundaries(
    tiny_tumor_bam: str,
) -> None:
    with pysam.AlignmentFile(tiny_tumor_bam) as bam:
        recording_bam = RecordingAlignmentFile(bam)
        candidates = [
            (("1", 8450), "T", "G"),
            (("1", 8450), "T", "A"),
            (("1", 8492), "G", "A"),
        ]

        batch_features = collect_bam_features_batch(
            cast(Any, recording_bam),
            candidates,
            max_cluster_span=40,
        )

        assert len(recording_bam.fetch_calls) == 2
        for candidate in candidates:
            expected = BamFeatures.from_alignment_file(
                bam,
                candidate[0],
                candidate[1],
                candidate[2],
            )
            assert_bam_features_equal(
                batch_features[(candidate[0], candidate[1], candidate[2])],
                expected,
            )
