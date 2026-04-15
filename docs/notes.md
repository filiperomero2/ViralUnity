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

2. **Reference extraction** — for each selected sample, a unique assembly key (`ref_key`) is constructed from the matched family and the reference accession (`{family}_{accession}`). The corresponding genome is extracted from `--viral-genomes`. The selection strategy controls how the reference accession is chosen (see below).

3. **Consensus assembly** — uses the same alignment + consensus rules as `viralunity consensus`, with the dynamically-selected reference. Output goes to `assembly/{ref_key}/`.

### Reference selection strategies

#### `--reference-selection-strategy taxid` (default)

The taxid assigned to each read or contig by the classifier (Kraken2 or Diamond) is looked up in `--viral-taxids` (`genome2taxid.tsv`) to find all genome accessions that share that taxid. Every matching accession generates a separate reference assembly for that sample.

**When to use:** when your classifier database and your `--viral-genomes` database were both built from the same RefSeq release. The taxid linkage is direct and requires no sequence comparison. Fast and deterministic.

**Limitation:** a taxid match does not guarantee the reference is the closest relative in the database — it only guarantees taxonomic identity. For highly diverse families (e.g. Flaviviridae), multiple accessions may match the same taxid, each generating its own assembly.

**Required databases:** `--viral-genomes` and `--viral-taxids` (both produced by `viralunity get-databases virus-genome`).

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-reference-assembly \
    --reference-selection-strategy taxid \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

#### `--reference-selection-strategy similarity`

De novo assembled contigs from each sample are BLASTed against `--viral-genomes` using `blastn`. The best hit that passes `--blast-qcov` and `--blast-pident` thresholds is selected as the reference.

**When to use:** when sequence identity to a specific genome matters more than taxonomic label — for example, when working with divergent or novel strains where taxid-based lookup may return a reference that is genetically distant. Requires de novo assembly (`--run-denovo-assembly`) to produce the contigs used as BLAST queries.

**Limitation:** requires a pre-built BLAST index alongside `--viral-genomes` (created automatically by `viralunity get-databases virus-genome`). If no contigs pass the identity/coverage thresholds, no reference assembly is triggered for that sample.

**Required databases:** `--viral-genomes` with its BLAST index (produced by `viralunity get-databases virus-genome`).

```bash
viralunity meta illumina \
    --sample-sheet samples.csv \
    --config-file config.yaml \
    --output /path/to/output \
    --kraken2-database /path/to/kraken2_db \
    --krona-database /path/to/krona_taxonomy \
    --taxdump /path/to/taxdump \
    --run-denovo-assembly \
    --run-reference-assembly \
    --reference-selection-strategy similarity \
    --method kraken2 \
    --source reads \
    --families Coronaviridae,Orthomyxoviridae \
    --reads-count 100 \
    --blast-qcov 80 \
    --blast-pident 80 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --threads 4 --threads-total 8
```

### Strategy comparison

| | `taxid` | `similarity` |
|---|---|---|
| Requires de novo assembly | No | Yes |
| Requires `--viral-taxids` | Yes | No |
| Requires BLAST index | No | Yes |
| Selection basis | Taxonomic ID match | Sequence identity |
| Speed | Fast | Slower (BLAST per sample) |
| Best for | Known strains with good RefSeq coverage | Divergent or novel strains |

### Required databases

| Option | File produced by | Used by |
|--------|-----------------|---------|
| `--viral-genomes` | `viralunity get-databases virus-genome` | Both strategies |
| `--viral-taxids` | `viralunity get-databases virus-genome` | `taxid` only |
| BLAST index (`.nhr`/`.nin`/`.nsq`) | `viralunity get-databases virus-genome` | `similarity` only |

## Running tests

```bash
pytest -v test
```
