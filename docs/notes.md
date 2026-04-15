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

## Dynamic reference assembly

When `--run-reference-assembly` is enabled, the metagenomics pipeline adds a post-classification step that automatically selects reference genomes and runs consensus assembly for each sample that meets the hit threshold.

### How it works

1. **Selection checkpoint** — after taxonomic classification, a Snakemake checkpoint reads the taxa summary TSVs and selects candidate references. The `--method` and `--source` flags control which classification results are used (e.g., `kraken2` on `reads`). Only samples with at least `--reads-count` reads (or `--contigs-count` contigs) assigned to one of the `--families` are selected.

2. **Reference extraction** — for each selected sample, a unique assembly key (`ref_key`) is constructed from the matched family and the reference accession (`{family}_{accession}`). The corresponding genome is extracted from `--viral-genomes`. Two selection strategies are available:
   - `taxid` (default): the taxid from the classification result is looked up in `--viral-taxids` to find the matching accession.
   - `similarity`: the de novo contigs are BLASTed against `--viral-genomes` and the best hit above `--blast-qcov` / `--blast-pident` is used.

3. **Consensus assembly** — uses the same alignment + consensus rules as `viralunity consensus`, with the dynamically-selected reference. Output goes to `assembly/{ref_key}/`.

### Required databases

| Option | File produced by |
|--------|-----------------|
| `--viral-genomes` | `viralunity get-databases virus-genome` |
| `--viral-taxids` | `viralunity get-databases virus-genome` |

### Example

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-reference-assembly \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

## Running tests

```bash
pytest -v test
```
