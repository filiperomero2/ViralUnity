#!/usr/bin/env python

import pandas as pd
import glob
import os
import sys

params = snakemake.params
output_file = snakemake.output[0]

# Load taxid-to-family mapping
taxid_to_family = pd.read_csv(params['taxid_to_family'], sep=',', dtype=str)
taxid_to_family = taxid_to_family.rename(columns={taxid_to_family.columns[0]: "taxid", taxid_to_family.columns[1]: "family"})

# Load all diamond krona input TSVs
#diamond_inputs = glob.glob(os.path.join(params['path'], "*.diamond.krona_input.tsv"))
diamond_inputs = glob.glob(os.path.join(params['path'], "*.diamond.supported.krona_input.tsv"))

summary_data = []

for path in diamond_inputs:
    sample = os.path.basename(path).replace(".diamond.krona_input.tsv", "")
    try:
        df = pd.read_csv(path, sep='\t', header=None, names=['read', 'taxid'], dtype=str)
        total_reads = len(df)

        if total_reads == 0:
            print(f"Warning: {path} is empty.", file=sys.stderr)
            continue  # Skip empty files entirely

        # Split into classified and unclassified
        classified = df[df['taxid'] != '0']
        unclassified_count = df['taxid'].value_counts().get('0', 0)

        family_counts_list = []

        if not classified.empty:
            # Merge classified with family names
            merged = classified.merge(taxid_to_family, on="taxid", how="left")
            family_counts = merged["family"].value_counts().reset_index()
            family_counts.columns = ["family", "count"]
            family_counts["sample"] = sample
            family_counts_list.append(family_counts)

        if unclassified_count > 0:
            unclassified_row = pd.DataFrame([{
                "family": "unclassified",
                "count": unclassified_count,
                "sample": sample
            }])
            family_counts_list.append(unclassified_row)

        if family_counts_list:
            family_counts = pd.concat(family_counts_list, ignore_index=True)
            # Calculate percentage
            family_counts["percent"] = (family_counts["count"] / total_reads * 100).round(2)
            summary_data.append(family_counts)

    except Exception as e:
        print(f"Error processing {path}: {e}", file=sys.stderr)


# Combine and write
if summary_data:
    summary_df = pd.concat(summary_data, ignore_index=True)
    
    # Compute per-family maximum count across all samples (excluding unclassified)
    max_counts = summary_df[summary_df["family"] != "unclassified"].groupby("family")["count"].max()

    # Annotate observations
    def annotate(row):
        if row["family"] == "unclassified":
            return ""
        max_count = max_counts.get(row["family"], 0)
        if max_count > 0 and row["count"] < max_count * 0.01:
            return f"Below 1% threshold for {row['family']}."
        return ""

    summary_df["observation"] = summary_df.apply(annotate, axis=1)

    # Sort for presentation
    summary_df = summary_df.sort_values(by=["sample", "percent"], ascending=[True, False])

    # Select columns
    summary_df = summary_df[['sample', 'family', 'count', 'percent', 'observation']]
    summary_df.to_csv(output_file, index=False)
else:
    pd.DataFrame(columns=["sample", "family", "count", "percent", "observation"]).to_csv(output_file, index=False)
