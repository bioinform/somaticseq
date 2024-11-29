import click

from somaticseq.genomic_file_parsers.concat import main as concat
from somaticseq.run_somaticseq import main as somaticseq_single_thread
from somaticseq.single_sample_vcf2tsv import main as single_sample_vcf2tsv
from somaticseq.somatic_tsv2vcf import main as somatic_tsv2vcf
from somaticseq.somatic_vcf2tsv import main as somatic_vcf2tsv
from somaticseq.somatic_xgboost import main as somatic_xgboost
from somaticseq.somaticseq_parallel import main as somaticseq_multi_thread
from somaticseq.utilities.dockered_pipelines.makeAlignmentScripts import (
    main as make_alignment_scripts,
)
from somaticseq.utilities.dockered_pipelines.makeSomaticScripts import (
    main as make_somatic_scripts,
)
from somaticseq.utilities.dockered_pipelines.run_workflows import main as run_workflows
from somaticseq.utilities.linguistic_sequence_complexity import (
    main as linguistic_sequence_complexity,
)
from somaticseq.utilities.lociCounterWithLabels import main as loci_counter_with_labels
from somaticseq.utilities.paired_end_bam2fastq import main as paired_end_bam2fastq
from somaticseq.utilities.split_bed_into_equal_regions import (
    main as split_bed_into_equal_regions,
)
from somaticseq.vcf_modifier.splitVcf import main as split_vcf


@click.group()
def main() -> None:
    """Consensam commands."""
