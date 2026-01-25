#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
from typing import List


def infer_group_cols(df: pd.DataFrame, extra_group_cols: List[str] | None = None) -> List[str]:
    """
    Decide what defines a 'taxon' for max-RPM purposes.

    Always include rank+taxid if present.
    Optionally include tool/mode if they exist in the table (keeps categories separate).
    Optionally include user-provided extra columns (e.g. database, classifier, etc.).
    """
    group_cols = []

    # Keep tool/mode separation if present (matches your earlier design: apply per category)
    for c in ["tool", "mode"]:
        if c in df.columns:
            group_cols.append(c)

    # Rank and taxid identify the taxon at a given level
    if "rank" in df.columns:
        group_cols.append("rank")
    if "taxid" in df.columns:
        group_cols.append("taxid")
    else:
        raise ValueError("Input must contain 'taxid' column to identify taxa.")

    if extra_group_cols:
        for c in extra_group_cols:
            if c not in df.columns:
                raise ValueError(f"Requested group column '{c}' not found in input columns.")
            group_cols.append(c)

    return group_cols


def apply_bleed_filter(
    df: pd.DataFrame,
    rpm_col: str = "rpm",
    fraction: float = 0.005,
    rpm_floor: float = 1.0,
    group_cols: List[str] | None = None,
) -> pd.DataFrame:
    """
    Add max-RPM bleed-through filter columns.

    - max_rpm: maximum RPM for the group across samples
    - bleed_threshold: fraction * max_rpm (only when max_rpm >= rpm_floor)
    - bleed_pass: True if rpm >= bleed_threshold OR max_rpm < rpm_floor (not applied)
    - bleed_applied: True/False indicating whether filter was applied for that taxon
    """
    if "sample" not in df.columns:
        raise ValueError("Input must contain 'sample' column.")
    if rpm_col not in df.columns:
        raise ValueError(f"Input missing RPM column '{rpm_col}'. Did you run add_RPM_to_summary first?")

    out = df.copy()
    if group_cols is None:
        group_cols = infer_group_cols(out)

    # Compute max rpm per taxon group
    max_rpm = out.groupby(group_cols, dropna=False)[rpm_col].max().reset_index().rename(columns={rpm_col: "max_rpm"})
    out = out.merge(max_rpm, on=group_cols, how="left")

    # Determine threshold and pass/fail
    out["bleed_applied"] = out["max_rpm"] >= float(rpm_floor)

    # Default threshold is 0 (so everything passes) when not applied; but keep explicit threshold column too
    out["bleed_threshold"] = 0.0
    mask = out["bleed_applied"]
    out.loc[mask, "bleed_threshold"] = out.loc[mask, "max_rpm"].astype(float) * float(fraction)

    # Pass if above threshold, OR if not applied (max_rpm below floor)
    out["bleed_pass"] = True
    out.loc[mask, "bleed_pass"] = out.loc[mask, rpm_col].astype(float) >= out.loc[mask, "bleed_threshold"].astype(float)

    return out


def run_cli():
    ap = argparse.ArgumentParser(
        description="Apply max-RPM bleed-through filter to an RPM-augmented taxa summary TSV."
    )
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (must include sample,taxid,rpm).")
    ap.add_argument("--out", required=True, help="Output TSV with bleed filter columns added.")
    ap.add_argument("--rpm-col", default="rpm", help="Name of RPM column (default: rpm).")
    ap.add_argument("--fraction", type=float, default=0.005, help="Threshold fraction of max_rpm (default: 0.005).")
    ap.add_argument(
        "--rpm-floor",
        type=float,
        default=1.0,
        help="If max_rpm < this floor, do NOT apply bleed filter for that taxon (default: 1.0).",
    )
    ap.add_argument(
        "--group-cols",
        default=None,
        help="Optional comma-separated columns to define groups. If omitted, uses tool/mode (if present) + rank + taxid.",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t")

    group_cols = None
    if args.group_cols:
        group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]

    out = apply_bleed_filter(
        df,
        rpm_col=args.rpm_col,
        fraction=args.fraction,
        rpm_floor=args.rpm_floor,
        group_cols=group_cols,
    )
    out.to_csv(args.out, sep="\t", index=False)


def run_snakemake():
    inp = str(snakemake.input[0])
    outp = str(snakemake.output[0])

    rpm_col = getattr(snakemake.params, "rpm_col", "rpm")
    fraction = float(getattr(snakemake.params, "fraction", 0.005))
    rpm_floor = float(getattr(snakemake.params, "rpm_floor", 1.0))

    group_cols = getattr(snakemake.params, "group_cols", None)
    if group_cols is not None:
        # allow passing list or comma-separated string
        if isinstance(group_cols, str):
            group_cols = [c.strip() for c in group_cols.split(",") if c.strip()]
        elif not isinstance(group_cols, list):
            raise ValueError("snakemake.params.group_cols must be a list or comma-separated string.")

    df = pd.read_csv(inp, sep="\t")
    out = apply_bleed_filter(
        df,
        rpm_col=rpm_col,
        fraction=fraction,
        rpm_floor=rpm_floor,
        group_cols=group_cols,
    )
    out.to_csv(outp, sep="\t", index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
