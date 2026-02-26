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
    TRIM = "trim"
    TRIM_HEAD = "trim_head"
    TRIM_TAIL = "trim_tail"
    CUT_FRONT_MEAN_QUALITY = "cut_front_mean_quality"
    CUT_TAIL_MEAN_QUALITY = "cut_tail_mean_quality"
    CUT_RIGHT_WINDOW_SIZE = "cut_right_window_size"
    CUT_RIGHT_MEAN_QUALITY = "cut_right_mean_quality"
    REMOVE_HUMAN_READS = "remove_human_reads"
    REMOVE_UNCLASSIFIED_READS = "remove_unclassified_reads"
    # Extended metagenomics
    HOST_REFERENCE = "host_reference"
    TAXDUMP = "taxdump"
    RUN_DENOVO_ASSEMBLY = "run_denovo_assembly"
    RUN_KRAKEN2_READS = "run_kraken2_reads"
    RUN_KRAKEN2_CONTIGS = "run_kraken2_contigs"
    RUN_DIAMOND_READS = "run_diamond_reads"
    RUN_DIAMOND_CONTIGS = "run_diamond_contigs"
    ASSEMBLY_SUMMARY = "assembly_summary"
    DIAMOND_DATABASE = "diamond_database"
    DIAMOND_SENSITIVITY = "diamond_sensitivity"
    EVALUE = "evalue"
    BLEED_FRACTION = "bleed_fraction"
    NEGATIVE_CONTROLS = "negative_controls"
    NEGATIVE_P_THRESHOLD = "negative_p_threshold"


class SampleSheetPattern:
    """Patterns for identifying sample files."""
    R1 = "R1"
    BARCODE = "barcode"


class SampleSheetSeparator:
    """Separators for parsing sample names."""
    UNDERSCORE = "_"
    HYPHEN = "-"
    DOT = "."

