"""Custom exceptions for ViralUnity pipeline."""


class ViralUnityError(Exception):
    """Base exception for all ViralUnity errors."""

    pass


class ValidationError(ViralUnityError):
    """Raised when validation of arguments or data fails."""

    pass


class ViralUnityFileNotFoundError(ViralUnityError):
    """Raised when a required file or directory is not found."""

    pass


class ConfigurationError(ViralUnityError):
    """Raised when there's an error in configuration."""

    pass


class SampleSheetError(ViralUnityError):
    """Raised when there's an error processing the sample sheet."""

    pass


class SampleConfigurationNotFoundError(ViralUnityError):
    """Raised when the sample configuration is not found."""

    pass


class Kraken2DatabaseNotFoundError(ViralUnityError):
    """Raised when the Kraken2 database is not found."""

    pass


class KronaDatabaseNotFoundError(ViralUnityError):
    """Raised when the Krona database is not found."""

    pass


class ReferenceNotFoundError(ViralUnityError):
    """Raised when the reference sequence file is not found."""

    pass


class PrimerSchemeNotFoundError(ViralUnityError):
    """Raised when the primer scheme file is not found."""

    pass


class AdaptersNotFoundError(ViralUnityError):
    """Raised when the Illumina adapter sequences file is not found or not provided."""

    pass


class TaxdumpNotFoundError(ViralUnityError):
    """Raised when the NCBI taxdump directory is not found or invalid."""

    pass


class DiamondDatabaseNotFoundError(ViralUnityError):
    """Raised when the Diamond database or assembly summary is not found."""

    pass
