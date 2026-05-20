"""Validation functions for ViralUnity pipeline arguments and data."""

import os
from typing import Dict, Any, List
import pandas as pd

from viralunity.exceptions import (
    Kraken2DatabaseNotFoundError,
    KronaDatabaseNotFoundError,
    ValidationError,
    FileNotFoundError,
    SampleSheetError,
    SampleConfigurationNotFoundError,
    ReferenceNotFoundError,
    PrimerSchemeNotFoundError,
    TaxdumpNotFoundError,
    DiamondDatabaseNotFoundError,
    AdaptersNotFoundError,
)
from viralunity.constants import DataType


def validate_file_exists(file_path: str, description: str = "File") -> None:
    """Validate that a file exists.

    Args:
        file_path: Path to the file
        description: Description of the file for error messages

    Raises:
        FileNotFoundError: If the file does not exist
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{description} does not exist: {file_path}")


def validate_directory_exists(
    directory_path: str, description: str = "Directory"
) -> None:
    """Validate that a directory exists.

    Args:
        directory_path: Path to the directory
        description: Description of the directory for error messages

    Raises:
        FileNotFoundError: If the directory does not exist
    """
    if not os.path.isdir(directory_path):
        raise FileNotFoundError(f"{description} does not exist: {directory_path}")


def validate_sample_sheet(
    sample_sheet_path: str, data_type: str
) -> Dict[str, List[str]]:
    """Validate and parse sample sheet file.

    Args:
        sample_sheet_path: Path to the sample sheet CSV file
        data_type: Type of sequencing data (illumina or nanopore)

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        SampleSheetError: If the sample sheet is invalid
        FileNotFoundError: If the sample sheet file does not exist
    """
    validate_file_exists(sample_sheet_path, "Sample sheet file")

    try:
        df = pd.read_csv(sample_sheet_path, header=None)
    except Exception as e:
        raise SampleSheetError(f"Failed to read sample sheet: {e}") from e

    samples = {}

    for idx in df.index:
        sample_name = df.iloc[idx, 0]

        if data_type == DataType.ILLUMINA:
            if df.shape[1] < 3:
                raise SampleSheetError(
                    f"Illumina sample sheet must have at least 3 columns. "
                    f"Found {df.shape[1]} columns."
                )

            sample_r1 = df.iloc[idx, 1]
            sample_r2 = df.iloc[idx, 2]

            validate_file_exists(sample_r1, f"R1 file for sample {sample_name}")
            validate_file_exists(sample_r2, f"R2 file for sample {sample_name}")

            samples[sample_name] = [sample_r1, sample_r2]

        else:  # nanopore
            if df.shape[1] < 2:
                raise SampleSheetError(
                    f"Nanopore sample sheet must have at least 2 columns. "
                    f"Found {df.shape[1]} columns."
                )

            sample_file = df.iloc[idx, 1]
            validate_file_exists(sample_file, f"File for sample {sample_name}")

            samples[sample_name] = [sample_file]

    if not samples:
        raise SampleSheetError("No valid samples found in sample sheet")

    return samples


def validate_illumina_requirements(args: Dict[str, Any]) -> None:
    """Validate Illumina-specific requirements (fastp: adapters optional)."""
    if args.get("data_type") != DataType.ILLUMINA:
        return

    adapters = args.get("adapters")
    if adapters and str(adapters).strip() and str(adapters).strip() != "NA":
        try:
            validate_file_exists(adapters, "Illumina adapter sequences file")
        except FileNotFoundError as e:
            raise AdaptersNotFoundError(
                f"Illumina adapter sequences file does not exist: {e}"
            ) from e


def validate_nanopore_requirements(args: Dict[str, Any]) -> None:
    """Validate Nanopore-specific requirements (e.g. polishing options)."""
    if args.get("data_type") != DataType.NANOPORE:
        return
    # Optional: validate medaka_model if provided (Medaka accepts known model names)
    # For now no strict validation; extend as needed.


def validate_consensus_requirements(args: Dict[str, Any]) -> None:
    """Validate consensus pipeline requirements.

    Exactly one of --reference or --segmented-reference must be provided.
    When --segmented-reference is used, the values are parsed from
    SEGMENT=PATH format and stored as a dict in args["reference"].

    Args:
        args: Dictionary of pipeline arguments

    Raises:
        ValidationError: If consensus requirements are not met
    """
    reference = args.get("reference")
    segmented_reference = args.get("segmented_reference")

    if reference and segmented_reference:
        raise ValidationError(
            "--reference and --segmented-reference are mutually exclusive. "
            "Please provide only one."
        )

    if not reference and not segmented_reference:
        raise ValidationError(
            "A reference is required. Provide --reference for a single reference "
            "or --segmented-reference for segmented viruses."
        )

    if segmented_reference:
        if isinstance(segmented_reference, dict):
            parsed_segments = segmented_reference
        else:
            parsed_segments = {}
            for entry in segmented_reference:
                if "=" not in entry:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Expected SEGMENT=PATH (e.g. S=/path/to/S.fasta)"
                    )
                segment_name, segment_path = entry.split("=", 1)
                segment_name = segment_name.strip()
                segment_path = segment_path.strip()
                if not segment_name or not segment_path:
                    raise ValidationError(
                        f"Invalid segmented reference format: '{entry}'. "
                        f"Both segment name and path are required."
                    )
                parsed_segments[segment_name] = segment_path

        args["reference"] = parsed_segments
        args["segmented_reference"] = None
        reference = parsed_segments

    if isinstance(reference, dict):
        for segment_name, segment_path in reference.items():
            try:
                validate_file_exists(
                    segment_path, f"Reference file for segment '{segment_name}'"
                )
            except FileNotFoundError as e:
                raise ReferenceNotFoundError(str(e)) from e
    else:
        try:
            validate_file_exists(reference, "Reference sequence file")
        except FileNotFoundError as e:
            raise ReferenceNotFoundError(f"Reference sequence file does not exist: {e}") from e

    primer_scheme = args.get("primer_scheme")
    if primer_scheme:
        try:
            validate_file_exists(primer_scheme, "Primer scheme file")
        except FileNotFoundError as e:
            raise PrimerSchemeNotFoundError(f"Primer scheme file does not exist: {e}") from e


def validate_metagenomics_requirements(args: Dict[str, Any]) -> None:
    """Validate metagenomics pipeline requirements.

    Kraken2 and Diamond are optional; validate only the resources for the tools
    the user has enabled (run_kraken2_reads, run_kraken2_contigs, run_diamond_reads,
    run_diamond_contigs).
    """
    run_denovo = args.get("run_denovo_assembly", True)
    run_k2_reads = args.get("run_kraken2_reads", True)
    run_k2_contigs = args.get("run_kraken2_contigs", True)
    run_diamond_reads = args.get("run_diamond_reads", False)
    run_diamond_contigs = args.get("run_diamond_contigs", False)
    
    if not run_denovo:
        if run_k2_contigs:
            raise ValidationError(
                "Cannot run kraken2 on contigs (--run-kraken2-contigs) when denovo assembly is disabled (--no-denovo-assembly)."
            )
        if run_diamond_contigs:
            raise ValidationError(
                "Cannot run diamond on contigs (--run-diamond-contigs) when denovo assembly is disabled (--no-denovo-assembly)."
            )

    any_kraken2 = run_k2_reads or run_k2_contigs
    any_diamond = run_diamond_reads or run_diamond_contigs
    any_classification = any_kraken2 or any_diamond

    # Krona database: required whenever any classification is run (Krona plots for both tools)
    if any_classification:
        krona_db = args.get("krona_database")
        if not krona_db or krona_db == "NA":
            raise KronaDatabaseNotFoundError(
                "Krona database directory is required when running any classification "
                "(Kraken2 and/or Diamond). Set --krona-database."
            )
        try:
            validate_directory_exists(krona_db, "Krona database directory")
        except FileNotFoundError as e:
            raise KronaDatabaseNotFoundError(
                f"Krona database directory does not exist: {e}"
            ) from e

    # Kraken2 database: required only when Kraken2 is enabled
    if any_kraken2:
        kraken2_db = args.get("kraken2_database")
        if not kraken2_db or kraken2_db == "NA":
            raise Kraken2DatabaseNotFoundError(
                "Kraken2 database directory is required when running Kraken2. "
                "Set --kraken2-database or disable Kraken2 with --no-kraken2-reads / --no-kraken2-contigs."
            )
        try:
            validate_directory_exists(kraken2_db, "Kraken2 database directory")
        except FileNotFoundError as e:
            raise Kraken2DatabaseNotFoundError(
                f"Kraken2 database directory does not exist: {e}"
            ) from e

    # Taxdump: required for taxonomic summaries whenever any classification is run
    if any_classification:
        taxdump = args.get("taxdump", "NA")
        if not taxdump or taxdump == "NA":
            raise TaxdumpNotFoundError(
                "taxdump directory (NCBI nodes.dmp, names.dmp) is required for taxonomic summaries "
                "when running any classification. Set --taxdump."
            )
        try:
            validate_directory_exists(taxdump, "Taxdump directory")
        except FileNotFoundError as e:
            raise TaxdumpNotFoundError(f"Taxdump directory does not exist: {e}") from e
        nodes = os.path.join(taxdump, "nodes.dmp")
        names = os.path.join(taxdump, "names.dmp")
        if not os.path.isfile(nodes) or not os.path.isfile(names):
            raise TaxdumpNotFoundError(
                f"Taxdump directory must contain nodes.dmp and names.dmp: {taxdump}"
            )

    # Diamond: require database and assembly summary only when Diamond is enabled
    if any_diamond:
        diamond_db = args.get("diamond_database", "NA")
        if not diamond_db or diamond_db == "NA":
            raise DiamondDatabaseNotFoundError(
                "diamond_database is required when running Diamond. "
                "Set --diamond-database or do not use --run-diamond-reads / --run-diamond-contigs."
            )
        taxids = args.get("taxids", "NA")
        if not taxids or taxids == "NA":
            raise DiamondDatabaseNotFoundError(
                "taxids mapping file is required when running Diamond. "
                "Set --taxids or do not use --run-diamond-reads / --run-diamond-contigs."
            )
        if not os.path.isfile(taxids) and not (
            taxids.endswith(".gz") and os.path.isfile(taxids)
        ):
            raise DiamondDatabaseNotFoundError(
                f"Taxid mapping file not found: {taxids}"
            )

    # Deacon index: when provided for host depletion, must exist
    deacon_idx = args.get("deacon_index", "NA")
    if deacon_idx and str(deacon_idx).strip() not in ("", "NA"):
        try:
            validate_file_exists(deacon_idx, "Deacon index file")
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Deacon index file does not exist: {e}") from e

    validate_reference_assembly_requirements(args)


def validate_reference_assembly_requirements(args: Dict[str, Any]) -> None:
    """Validate cross-dependencies for Reference Assembly in metagenomics workflows."""
    if not args.get("run_reference_assembly"):
        return

    method = args.get("method")
    source = args.get("source")
    strategy = args.get("reference_selection_strategy")

    # If --similarity is used, it must be on contigs, and the respective classification tool must run on contigs
    if strategy == "similarity":
        run_denovo = args.get("run_denovo_assembly", False)
        if not run_denovo:
            raise ValidationError(
                "Reference selection strategy 'similarity' requires --run-denovo-assembly."
            )
        if source == "reads":
            raise ValidationError(
                "Reference selection strategy 'similarity' can only be used if --source includes 'contigs'."
            )

        # Check if the chosen method actually runs on contigs
        if method == "kraken2" and not args.get("run_kraken2_contigs", True):
            raise ValidationError(
                "Strategy 'similarity' with method 'kraken2' requires --run-kraken2-contigs."
            )
        if method == "diamond" and not args.get("run_diamond_contigs", False):
            raise ValidationError(
                "Strategy 'similarity' with method 'diamond' requires --run-diamond-contigs."
            )
        if method == "both" and not (
            args.get("run_kraken2_contigs", True)
            or args.get("run_diamond_contigs", False)
        ):
            raise ValidationError(
                "Strategy 'similarity' with method 'both' requires at least one mapping tool on contigs."
            )

    if method in ["kraken2", "both"]:
        if source in ["reads", "both"] and not args.get("run_kraken2_reads", True):
            raise ValidationError(
                "Method includes 'kraken2' on 'reads' but --no-kraken2-reads was passed."
            )
        if source in ["contigs", "both"] and not args.get("run_kraken2_contigs", True):
            raise ValidationError(
                "Method includes 'kraken2' on 'contigs' but --no-kraken2-contigs was passed."
            )

    if method in ["diamond", "both"]:
        if source in ["reads", "both"] and not args.get("run_diamond_reads", False):
            raise ValidationError(
                "Method includes 'diamond' on 'reads' but --run-diamond-reads is not enabled."
            )
        if source in ["contigs", "both"] and not args.get("run_diamond_contigs", False):
            raise ValidationError(
                "Method includes 'diamond' on 'contigs' but --run-diamond-contigs is not enabled."
            )

    viral_genomes = args.get("viral_genomes")
    if viral_genomes and not os.path.isfile(viral_genomes):
        raise FileNotFoundError(f"Viral genomes file does not exist: {viral_genomes}")

    if strategy == "taxid":
        viral_taxids = args.get("viral_taxids")
        if viral_taxids and not os.path.isfile(viral_taxids):
            raise FileNotFoundError(f"Viral taxids file does not exist: {viral_taxids}")


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
#
# CLI arguments that point to a filesystem location. The lists are kept here
# (next to the validators) so they can be reused both at validation time and
# inside ``viralunity_meta.main`` / ``viralunity_consensus.main`` to make the
# paths absolute before anything else sees them. Resolving paths up-front
# means a user who runs ``viralunity meta nanopore ... --host-reference
# databases/host/host.fasta --config-file scratch/run.yml`` does NOT silently
# get the host reference looked up under ``scratch/databases/host/...`` just
# because the generated config happens to live in ``scratch/``.

META_PATH_ARG_KEYS = (
    "sample_sheet",
    "config_file",
    "output",
    "kraken2_database",
    "krona_database",
    "taxdump",
    "host_reference",
    "deacon_index",
    "taxids",
    "diamond_database",
    "viral_genomes",
    "viral_taxids",
    "adapters",
)

CONSENSUS_PATH_ARG_KEYS = (
    "sample_sheet",
    "config_file",
    "output",
    "reference",
    "primer_scheme",
    "adapters",
)


def _is_path_sentinel(value: Any) -> bool:
    """Return True if a path-typed argument value should be left untouched.

    Sentinels are: ``None``, non-string scalars (ints, bools), the empty
    string, and the literal placeholder ``"NA"`` (after stripping
    whitespace). This matches the conventions used by the existing
    validators in this module.
    """
    if value is None:
        return True
    if not isinstance(value, str):
        return True
    stripped = value.strip()
    if not stripped:
        return True
    if stripped == "NA":
        return True
    return False


def resolve_path_args(
    args: Dict[str, Any],
    keys,
    base_dir: str = None,
) -> Dict[str, Any]:
    """Rewrite path-typed argument values to absolute paths in place.

    Each value listed in ``keys`` that is a non-sentinel relative path string
    is replaced with ``os.path.abspath(os.path.join(base_dir, value))``.
    Absolute paths, sentinels (``None``, ``""``, ``"NA"``), and non-string
    values are left unchanged. Missing keys are ignored.

    The ``reference`` argument of the consensus pipeline can be a dict
    (segmented reference, ``segment -> path``); each value of the dict is
    resolved while preserving the keys.

    Args:
        args: Mutable dict of CLI arguments.
        keys: Iterable of argument keys whose values are filesystem paths.
        base_dir: Base directory to resolve relative paths against.
            Defaults to the current working directory at call time.

    Returns:
        The same ``args`` dict (modified in place). Returned for chaining.
    """
    base_dir = base_dir or os.getcwd()

    for key in keys:
        if key not in args:
            continue
        value = args[key]

        if isinstance(value, dict):
            args[key] = {
                seg: (
                    v
                    if _is_path_sentinel(v) or os.path.isabs(v)
                    else os.path.abspath(os.path.join(base_dir, v))
                )
                for seg, v in value.items()
            }
            continue

        if _is_path_sentinel(value):
            continue
        if os.path.isabs(value):
            continue

        args[key] = os.path.abspath(os.path.join(base_dir, value))

    return args


def get_samples_from_args(args: Dict[str, Any]) -> Dict[str, List[str]]:
    """Extract and validate samples from arguments.

    Args:
        args: Dictionary of pipeline arguments

    Returns:
        Dictionary mapping sample names to file paths

    Raises:
        ValidationError: If samples cannot be determined from arguments
    """
    sample_sheet = args.get("sample_sheet")
    samples = args.get("samples")
    data_type = args.get("data_type")

    if sample_sheet and os.path.isfile(sample_sheet):
        return validate_sample_sheet(sample_sheet, data_type)
    elif samples:
        return samples
    else:
        raise SampleConfigurationNotFoundError(
            "Either 'sample_sheet' or 'samples' must be provided. "
            f"Sample sheet path: {sample_sheet}"
        )
