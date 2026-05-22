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
    REMOVE_HUMAN_READS = "remove_human_reads"
    REMOVE_UNCLASSIFIED_READS = "remove_unclassified_reads"
    # Extended metagenomics
    HOST_REFERENCE = "host_reference"
    DEACON_INDEX = "deacon_index"
    TAXDUMP = "taxdump"
    RUN_DENOVO_ASSEMBLY = "run_denovo_assembly"
    RUN_KRAKEN2_READS = "run_kraken2_reads"
    RUN_KRAKEN2_CONTIGS = "run_kraken2_contigs"
    RUN_DIAMOND_READS = "run_diamond_reads"
    RUN_DIAMOND_CONTIGS = "run_diamond_contigs"
    TAXIDS = "taxids"
    VIRAL_GENOMES = "viral_genomes"
    VIRAL_TAXIDS = "viral_taxids"
    DIAMOND_DATABASE = "diamond_database"
    DIAMOND_SENSITIVITY = "diamond_sensitivity"
    EVALUE = "evalue"
    MINIMUM_HIT_GROUP = "minimum_hit_group"
    DIAMOND_MAX_TARGET_SEQS = "diamond_max_target_seqs"
    KRAKEN2_EXTRA_FLAGS = "kraken2_extra_flags"
    MINIMAP2_CONSENSUS_ALIGN_FLAGS = "minimap2_consensus_align_flags"
    BLEED_FRACTION = "bleed_fraction"
    NEGATIVE_CONTROLS = "negative_controls"
    NEGATIVE_P_THRESHOLD = "negative_p_threshold"
    # Nanopore polishing
    RUN_POLISH_RACON = "run_polish_racon"
    RUN_POLISH_MEDAKA = "run_polish_medaka"
    MEDAKA_MODEL = "medaka_model"
    # Reference Assembly
    RUN_REFERENCE_ASSEMBLY = "run_reference_assembly"
    REF_ASSEMBLY_METHOD = "ref_assembly_method"
    REF_ASSEMBLY_SOURCE = "ref_assembly_source"
    REF_ASSEMBLY_READS_COUNT = "ref_assembly_reads_count"
    REF_ASSEMBLY_CONTIGS_COUNT = "ref_assembly_contigs_count"
    REF_ASSEMBLY_FAMILIES = "ref_assembly_families"
    REF_SELECTION_STRATEGY = "ref_selection_strategy"
    REF_BLAST_QCOV = "ref_blast_qcov"
    REF_BLAST_PIDENT = "ref_blast_pident"

    AF_THRESHOLD = "af_threshold"
    AF_ISNV_THRESHOLD = "af_isnv_threshold"
    REMOVE_HUMAN_READS = "remove_human_reads"
    REMOVE_UNCLASSIFIED_READS = "remove_unclassified_reads"
    CHUNK_SIZE = "chunk_size"
    CLAIR3_MODEL = "clair3_model"
    VARIANT_QUALITY = "variant_quality"
    VARIANT_DEPTH = "variant_depth"
    MINIMUM_MAP_QUALITY = "minimum_map_quality"
    RUN_ISNV = "run_isnv"


class SampleSheetPattern:
    """Patterns for identifying sample files."""

    R1 = "R1"
    BARCODE = "barcode"


class SampleSheetSeparator:
    """Separators for parsing sample names."""

    UNDERSCORE = "_"
    HYPHEN = "-"
    DOT = "."


class ResourceDefaults:
    """Default resource allocation for computational Snakemake rules.

    Each rule gets <rule>_cpus and <rule>_ram (GB) in the config YAML.
    """

    DEFAULT_CPUS = 2
    DEFAULT_RAM = 4  # GB

    # Consensus pipeline — Illumina computational rules
    CONSENSUS_ILLUMINA_RULES = [
        "perform_qc",
        "map_reads",
        "trim_primer_sequences",
        "detect_isnv",
    ]

    # Consensus pipeline — Nanopore computational rules
    CONSENSUS_NANOPORE_RULES = [
        "map_reads",
        "trim_primer_sequences",
        "infer_consensus_sequence",
    ]

    # Metagenomics pipeline — shared computational rules
    META_SHARED_RULES = [
        "remove_host_reads",
        "index_host_genome",
        "run_megahit",
        "run_kraken2_reads",
        "run_kraken2_contigs",
        "run_diamond_reads",
        "run_diamond_contigs",
        "index_viral_contigs",
        "remap_reads_to_viral_contigs",
        "bam_sort_index_idxstats_from_medaka",
    ]

    # Metagenomics pipeline — Illumina-specific computational rules
    META_ILLUMINA_RULES = [
        "perform_qc",
    ]

    # Metagenomics pipeline — Nanopore-specific computational rules
    META_NANOPORE_RULES = [
        "run_racon",
        "run_medaka",
    ]
