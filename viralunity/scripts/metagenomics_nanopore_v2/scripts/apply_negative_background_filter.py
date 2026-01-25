#!/usr/bin/env python3
import argparse
import math
from typing import List, Optional

import pandas as pd


def poisson_sf(k: int, mu: float) -> float:
    """Survival function P(X >= k) for X ~ Poisson(mu)."""
    if k <= 0:
        return 1.0
    if mu <= 0.0:
        return 0.0

    p = math.exp(-mu)  # pmf(0)
    cdf = p
    for i in range(1, k):
        p *= mu / i
        cdf += p
        if cdf > 1.0 - 1e-15:
            return 0.0

    sf = 1.0 - cdf
    return max(0.0, sf)


def infer_group_cols(df: pd.DataFrame, extra_group_cols: Optional[List[str]] = None) -> List[str]:
    """Default grouping: (tool,mode if present) + rank + taxid"""
    group_cols = []
    for c in ["tool", "mode"]:
        if c in df.columns:
            group_cols.append(c)

    if "rank" in df.columns:
        group_cols.append("rank")

    if "taxid" not in df.columns:
        raise ValueError("Input must contain 'taxid' column.")
    group_cols.append("taxid")

    if extra_group_cols:
        for c in extra_group_cols:
            if c not in df.columns:
                raise ValueError(f"Requested extra group column '{c}' not found.")
            group_cols.append(c)

    return group_cols


def apply_negative_background_filter(
    df: pd.DataFrame,
    negatives: List[str],
    count_col: str,
    total_reads_col: str = "total_reads",
    p_threshold: float = 0.01,
    group_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Negative-control background filter WITHOUT pseudocounts.

    For taxa present in negatives:
      b_t = sum(count_col) over negatives for taxon t
      B   = sum(total_reads) across negative samples (within category)
      lambda_t = b_t / B

      mu_{s,t} = lambda_t * total_reads_s
      p_bg = P(Poisson(mu_{s,t}) >= observed)

    For taxa absent in negatives:
      neg_evaluated=False, p_bg/mu_bg/neg_pass = NA
    """
    if "sample" not in df.columns:
        raise ValueError("Input must contain a 'sample' column.")
    if count_col not in df.columns:
        raise ValueError(f"Input missing count column '{count_col}'.")
    if total_reads_col not in df.columns:
        raise ValueError(
            f"Input missing '{total_reads_col}'. Run the RPM-adding script first so total_reads exists."
        )

    out = df.copy()

    if group_cols is None:
        group_cols = infer_group_cols(out)

    # If no negative controls were provided, do not apply this filter.
    # Return table with expected columns present (NA) so downstream steps remain stable.
    if negatives is None or len(negatives) == 0:
        out = df.copy()
        out["is_negative_control"] = False
        out["neg_controls_used"] = 0
        out["p_threshold"] = float(p_threshold)
        out["neg_b_t"] = pd.NA
        out["neg_B"] = pd.NA
        out["lambda_bg"] = pd.NA
        out["neg_evaluated"] = False
        out["mu_bg"] = pd.NA
        out["p_bg"] = pd.NA
        out["neg_pass"] = pd.NA
        return out

    neg_set = set(negatives)
    out["is_negative_control"] = out["sample"].isin(neg_set)

    neg_samples_present = sorted(set(out.loc[out["is_negative_control"], "sample"].unique()))
    if len(neg_samples_present) == 0:
        raise ValueError("None of the provided negative controls are present in the input table.")

    # Category columns = group_cols excluding rank+taxid (for per-tool/mode totals)
    category_cols = [c for c in group_cols if c not in ["rank", "taxid"]]
    if len(category_cols) == 0:
        out["_cat"] = "__all__"
        category_cols = ["_cat"]

    # neg_B per category: sum unique negative sample total_reads (avoid duplication across taxa rows)
    neg_sample_totals = (
        out.loc[out["is_negative_control"], category_cols + ["sample", total_reads_col]]
        .drop_duplicates(subset=category_cols + ["sample"])
    )
    B_by_cat = (
        neg_sample_totals.groupby(category_cols, dropna=False)[total_reads_col]
        .sum()
        .reset_index()
        .rename(columns={total_reads_col: "neg_B"})
    )

    # b_t per taxon group: sum counts across negatives
    b_by_tax = (
        out.loc[out["is_negative_control"]]
        .groupby(group_cols, dropna=False)[count_col]
        .sum()
        .reset_index()
        .rename(columns={count_col: "neg_b_t"})
    )

    # attach neg_B to each taxon group via category
    b_by_tax = b_by_tax.merge(B_by_cat, on=category_cols, how="left")

    # lambda only defined for taxa seen in negatives (b_by_tax rows)
    b_by_tax["lambda_bg"] = b_by_tax["neg_b_t"].astype(float) / b_by_tax["neg_B"].astype(float)

    # merge taxon-level background params to all rows (taxa absent in negatives -> NaN)
    out = out.merge(
        b_by_tax[group_cols + ["neg_b_t", "neg_B", "lambda_bg"]],
        on=group_cols,
        how="left",
    )

    # also merge neg_B by category so we still report it even when taxon absent
    out = out.merge(B_by_cat, on=category_cols, how="left", suffixes=("", "_cat"))
    out["neg_B"] = out["neg_B"].fillna(out["neg_B_cat"])
    out = out.drop(columns=["neg_B_cat"])

    # Evaluation flag: only evaluate taxa that exist in b_by_tax (lambda_bg not null)
    out["neg_evaluated"] = out["lambda_bg"].notna()

    # mu and p only for evaluated rows; else NA
    out["mu_bg"] = pd.NA
    out["p_bg"] = pd.NA
    out["neg_pass"] = pd.NA

    eval_mask = out["neg_evaluated"]

    # observed counts (integer)
    obs = out.loc[eval_mask, count_col].fillna(0).astype(float).round().astype(int)
    mu = (out.loc[eval_mask, "lambda_bg"].astype(float) * out.loc[eval_mask, total_reads_col].astype(float))

    out.loc[eval_mask, "mu_bg"] = mu

    pvals = [poisson_sf(int(k), float(m)) for k, m in zip(obs.tolist(), mu.tolist())]
    out.loc[eval_mask, "p_bg"] = pvals
    out.loc[eval_mask, "neg_pass"] = (pd.Series(pvals, index=out.index[eval_mask]) < float(p_threshold))

    out["neg_controls_used"] = len(neg_samples_present)
    out["p_threshold"] = float(p_threshold)

    # cleanup helper
    if "_cat" in out.columns:
        out = out.drop(columns=["_cat"])

    return out


def run_cli():
    ap = argparse.ArgumentParser(
        description="Apply negative-control Poisson background filter (skip taxa absent from negatives)."
    )
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (must include sample,taxid,total_reads).")
    ap.add_argument("--out", required=True, help="Output TSV with background columns added.")
    ap.add_argument("--negatives", required=True, help="Comma-separated negative control sample IDs.")
    ap.add_argument("--count-col", required=True, help="Numerator column (e.g. count or mapped_reads).")
    ap.add_argument("--total-reads-col", default="total_reads", help="Denominator column (default: total_reads).")
    ap.add_argument("--p-threshold", type=float, default=0.01, help="Pass if p_bg < this value (default: 0.01).")
    ap.add_argument(
        "--group-cols",
        default=None,
        help="Optional comma-separated group columns. Default: tool/mode if present + rank + taxid.",
    )
    args = ap.parse_args()

    negatives = [x.strip() for x in args.negatives.split(",") if x.strip()]
    group_cols = None
    if args.group_cols:
        group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]

    df = pd.read_csv(args.inp, sep="\t")
    out = apply_negative_background_filter(
        df,
        negatives=negatives,
        count_col=args.count_col,
        total_reads_col=args.total_reads_col,
        p_threshold=args.p_threshold,
        group_cols=group_cols,
    )
    out.to_csv(args.out, sep="\t", index=False, na_rep="NA")


def run_snakemake():
    inp = str(snakemake.input[0])
    outp = str(snakemake.output[0])

    negatives = getattr(snakemake.params, "negatives", None)
    if negatives is None:
        raise ValueError("snakemake.params must include 'negatives' (list or comma-separated string).")
    if isinstance(negatives, str):
        negatives = [x.strip() for x in negatives.split(",") if x.strip()]
    if not isinstance(negatives, list):
        raise ValueError("snakemake.params.negatives must be a list or comma-separated string.")

    count_col = getattr(snakemake.params, "count_col", None)
    if count_col is None:
        raise ValueError("snakemake.params must include 'count_col' (e.g., 'count' or 'mapped_reads').")

    total_reads_col = getattr(snakemake.params, "total_reads_col", "total_reads")
    p_threshold = float(getattr(snakemake.params, "p_threshold", 0.01))

    group_cols = getattr(snakemake.params, "group_cols", None)
    if group_cols is not None:
        if isinstance(group_cols, str):
            group_cols = [c.strip() for c in group_cols.split(",") if c.strip()]
        elif not isinstance(group_cols, list):
            raise ValueError("snakemake.params.group_cols must be a list or comma-separated string.")

    df = pd.read_csv(inp, sep="\t")
    out = apply_negative_background_filter(
        df,
        negatives=negatives,
        count_col=count_col,
        total_reads_col=total_reads_col,
        p_threshold=p_threshold,
        group_cols=group_cols,
    )
    out.to_csv(outp, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        run_cli()
