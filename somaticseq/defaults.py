from typing import Literal

MIN_MAPPING_QUALITY: int = 1
MIN_BASE_QUALITY: int = 5
MIN_CALLER: float = 0.5
PASS_SCORE: float = 0.5
LOWQUAL_SCORE: float = 0.1
HOMOZYGOUS_FRAC: float = 0.85
HETEROZYGOUS_FRAC: float = 0.01

SNV_TSV_SUFFIX: str = "sSNV.tsv"
INDEL_TSV_SUFFIX: str = "sINDEL.tsv"
SNV_VCF_SUFFIX: str = "sSNV.vcf"
INDEL_VCF_SUFFIX: str = "sINDEL.vcf"
ENSEMBLE_PREFIX: str = "Ensemble."
CONSENSUS_PREFIX: str = "Consensus."
CLASSIFIED_PREFIX: str = "SSeq.Classified."
TUMOR_NAME: str = "TUMOR"
NORMAL_NAME: str = "NORMAL"

ALGORITHM: Literal["xgboost", "ada"] = "xgboost"
DEFAULT_XGB_BOOST_ROUNDS: int = 500
DEFAULT_NUM_TREES_PREDICT: int = 100
