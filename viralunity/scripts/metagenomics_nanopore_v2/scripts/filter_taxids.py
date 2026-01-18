import pandas as pd
import os
import sys

input_file = snakemake.input[0]
output_file = snakemake.output[0]

exclude_taxids = set(map(str, snakemake.params.exclude_taxids))
taxid_column = snakemake.params.taxid_column
keep_columns = snakemake.params.get("keep_columns", None)

# ----------------------------
# Handle empty input gracefully
# ----------------------------
if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
    # Create empty output and exit cleanly
    open(output_file, "w").close()
    print(
        f"[filter_taxids] Empty input file: {input_file}. "
        f"Created empty output: {output_file}",
        file=sys.stderr,
    )
    sys.exit(0)

# ----------------------------
# Read TSV
# ----------------------------
try:
    df = pd.read_csv(input_file, sep="\t", header=None, dtype=str)
except pd.errors.EmptyDataError:
    open(output_file, "w").close()
    print(
        f"[filter_taxids] No data in input file: {input_file}. "
        f"Created empty output.",
        file=sys.stderr,
    )
    sys.exit(0)

# ----------------------------
# Select columns (e.g. Kraken2 cut -f)
# ----------------------------
if keep_columns is not None:
    df = df.iloc[:, keep_columns]

# If after slicing nothing remains, exit cleanly
if df.empty:
    open(output_file, "w").close()
    print(
        f"[filter_taxids] No rows after column selection. "
        f"Created empty output.",
        file=sys.stderr,
    )
    sys.exit(0)

# ----------------------------
# Resolve taxid column
# ----------------------------
if isinstance(taxid_column, int):
    taxid_colname = df.columns[taxid_column]
else:
    taxid_colname = taxid_column
    if taxid_colname not in df.columns:
        open(output_file, "w").close()
        print(
            f"[filter_taxids] TaxID column '{taxid_colname}' not found. "
            f"Created empty output.",
            file=sys.stderr,
        )
        sys.exit(0)

# ----------------------------
# Filter
# ----------------------------
if exclude_taxids:
    df = df[~df[taxid_colname].isin(exclude_taxids)]

# ----------------------------
# Write output (Krona expects no header)
# ----------------------------
df.to_csv(output_file, sep="\t", index=False, header=False)

print(
    f"[filter_taxids] Wrote {len(df)} rows to {output_file}",
    file=sys.stderr,
)
