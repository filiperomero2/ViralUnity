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
    REMOVE_HUMAN_READS = "remove_human_reads"
    REMOVE_UNCLASSIFIED_READS = "remove_unclassified_reads"


class SampleSheetPattern:
    """Patterns for identifying sample files."""
    R1 = "R1"
    BARCODE = "barcode"


class SampleSheetSeparator:
    """Separators for parsing sample names."""
    UNDERSCORE = "_"
    HYPHEN = "-"
    DOT = "."

