#!/usr/bin/env python3
"""Filter and annotate a DIAMOND outfmt-6 TSV with read-support data from
samtools idxstats.

- Merge read counts from idxstats to base contig IDs (handles sub-regions).
- Filter DIAMOND rows to keep only contigs with >= min_mapped reads.
- Annotate each retained row with total mapped_reads (last column).
"""

import argparse
from pathlib import Path
from typing import Dict, Set, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter + annotate DIAMOND TSVs using read counts from idxstats."
    )
    parser.add_argument("--idxstats", required=True, help="samtools idxstats TSV")
    parser.add_argument("--diamond", required=True, help="DIAMOND TSV to filter")
    parser.add_argument("--output", required=True, help="Output filtered TSV path")
    parser.add_argument(
        "--min-mapped",
        type=int,
        default=1,
        help="Minimum mapped reads required per contig (default: 1)",
    )
    return parser.parse_args()


def load_mapped_counts(idxstats_path: Path) -> Dict[str, int]:
    """Return {base_contig: total_mapped_reads} summed across sub-regions."""
    counts: Dict[str, int] = {}
    with idxstats_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            ref = parts[0]
            if ref == "*":
                continue
            try:
                mapped = int(parts[2])
            except (ValueError, IndexError):
                continue
            base = ref.split("|", 1)[0]
            counts[base] = counts.get(base, 0) + mapped
    return counts


def diamond_base_ids(diamond_path: Path) -> Set[str]:
    """Return set of base contig IDs in the DIAMOND TSV."""
    bases: Set[str] = set()
    with diamond_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            qseqid = line.split("\t", 1)[0]
            bases.add(qseqid.split("|", 1)[0])
    return bases


def filter_and_annotate(
    diamond_path: Path,
    output_path: Path,
    mapped_counts: Dict[str, int],
    min_mapped: int,
) -> Tuple[int, int]:
    """Write filtered & annotated DIAMOND TSV. Returns (rows_kept, unique_contigs_kept)."""
    kept_rows = 0
    kept_bases: Set[str] = set()
    header_written = False

    with diamond_path.open() as inp, output_path.open("w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                if not header_written:
                    out.write("# added_column: mapped_reads\n")
                    header_written = True
                continue
            if not line.strip():
                continue
            qseqid = line.split("\t", 1)[0]
            base = qseqid.split("|", 1)[0]
            mapped = mapped_counts.get(base, 0)
            if mapped < min_mapped:
                continue
            out.write(line.rstrip("\n") + f"\t{mapped}\n")
            kept_rows += 1
            kept_bases.add(base)

    return kept_rows, len(kept_bases)


def main() -> None:
    if "snakemake" in globals():
        idxstats_path = Path(snakemake.input.idxstats)
        diamond_path = Path(snakemake.input.diamond)
        output_path = Path(snakemake.output.filtered)
        min_mapped = int(getattr(snakemake.params, "min_mapped", 1))
    else:
        args = parse_args()
        idxstats_path = Path(args.idxstats)
        diamond_path = Path(args.diamond)
        output_path = Path(args.output)
        min_mapped = args.min_mapped

    if idxstats_path.exists() and idxstats_path.stat().st_size == 0:
        print("idxstats file is empty, creating empty diamond output.")
        output_path.touch()
        return

    mapped_counts = load_mapped_counts(idxstats_path)
    total_bases = len(diamond_base_ids(diamond_path))

    kept_rows, kept_contigs = filter_and_annotate(
        diamond_path, output_path, mapped_counts, min_mapped
    )

    print(
        "[✓] DIAMOND filter/annotate summary:\n"
        f"    • unique contigs in original : {total_bases}\n"
        f"    • contigs retained          : {kept_contigs}\n"
        f"    • DIAMOND rows retained     : {kept_rows}\n"
        f"    → {output_path}\n"
    )


if __name__ == "__main__":
    main()
