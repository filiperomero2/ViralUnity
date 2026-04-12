# Notes

## Segmented viruses

ViralUnity natively supports the assembly of segmented viral genomes. Instead of running the pipeline multiple times, pass multiple segment references with `--segmented-reference`:

```bash
viralunity consensus illumina \
    --sample-sheet samples.csv \
    --config-file config_segmented.yml \
    --output /path/to/output \
    --segmented-reference S=/path/to/s_segment.fasta \
    --segmented-reference L=/path/to/l_segment.fasta
```

When this argument is detected, the pipeline dispatches a specialised modular workflow that processes each segment independently in parallel. Results are organised under `samples/{sample_name}/{segment_name}/`.

## Nanopore data

The metagenomics pipeline supports Nanopore data via `viralunity meta nanopore`. There is no fastp-based QC step; optional dehosting, MEGAHIT assembly, and Racon/Medaka polishing are available. Sample sheets have two columns (sample ID, path to one FASTQ/FASTA per sample).

## Reference header sanitization

The nanopore consensus workflow automatically sanitizes reference FASTA headers before use. Special characters (`/`, `\`, `|`, `,`, `~`, and spaces) in sequence identifiers are replaced with underscores (`_`). This prevents issues with downstream tools such as clair3 that use the sequence ID to create output directories.

## Running tests

```bash
pytest -v test
```
