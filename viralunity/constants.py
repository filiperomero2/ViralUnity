"""Constants used throughout the ViralUnity pipeline."""


class DataType:
    """Sequencing data types supported by ViralUnity."""
    ILLUMINA = "illumina"
    NANOPORE = "nanopore"
    
    @classmethod
    def is_valid(cls, value: str) -> bool:
        """Check if a value is a valid data type."""
        return value in [cls.ILLUMINA, cls.NANOPORE]


class ConfigKeys:
    """Keys used in configuration files."""
    SAMPLES = "samples"
    DATA = "data"
    OUTPUT = "output"
    THREADS = "threads"
    KRAKEN2_DATABASE = "kraken2_database"
    KRONA_DATABASE = "krona_database"
    REFERENCE = "reference"
    SCHEME = "scheme"
    ADAPTERS = "adapters"
    MINIMUM_LENGTH = "minimum_length"
    MINIMUM_DEPTH = "minimum_depth"
    TRIM_HEAD = "trim_head"
    TRIM_TAIL = "trim_tail"
    CUT_FRONT_MEAN_QUALITY = "cut_front_mean_quality"
    CUT_TAIL_MEAN_QUALITY = "cut_tail_mean_quality"
    CUT_RIGHT_WINDOW_SIZE = "cut_right_window_size"
    CUT_RIGHT_MEAN_QUALITY = "cut_right_mean_quality"
    AF_THRESHOLD = "af_threshold"
    AF_ISNV_THRESHOLD = "af_isnv_threshold"
    REMOVE_HUMAN_READS = "remove_human_reads"
    REMOVE_UNCLASSIFIED_READS = "remove_unclassified_reads"
    CHUNK_SIZE = "chunk_size"
    CLAIR3_MODEL = "clair3_model"
    VARIANT_QUALITY = "variant_quality"
    MINIMUM_MAP_QUALITY = "minimum_map_quality"

class SampleSheetPattern:
    """Patterns for identifying sample files."""
    R1 = "R1"
    BARCODE = "barcode"


class SampleSheetSeparator:
    """Separators for parsing sample names."""
    UNDERSCORE = "_"
    HYPHEN = "-"
    DOT = "."

