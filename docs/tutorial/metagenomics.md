# Metagenomics pipeline tutorial

The metagenomics pipeline (`viralunity meta`) classifies reads of unknown content, optionally assembles them, and — if you ask it to — automatically picks reference genomes from the classification hits and runs a per-virus consensus assembly on the fly. This page builds it up incrementally on the same SARS-CoV-2 data you used in [the consensus tutorial](consensus.md), starting with a minimal Kraken2 run and ending with the full pipeline plus reference assembly.

```{note}
This page assumes you have completed [Setup](setup.md), including downloading the metagenomic databases under `databases/`.
```

## When to use it

Use `viralunity meta` when the contents of your sample are unknown, when you want to screen for unexpected viruses, when you need a community profile across many samples, or when you want both a classification and per-virus consensus genomes in one workflow. If you already know what is in the sample, the [consensus pipeline](consensus.md) is simpler and faster.

## What it does

The full pipeline strings together up to five conceptual phases, most of them optional and toggled by CLI flags. Anything in `[brackets]` is an optional step that only runs when you ask for it.

```text
raw FASTQ
  ──► [fastp QC, Illumina only]
  ──► [host depletion: Deacon  OR  minimap2 vs --host-reference]
  ├─► Kraken2 on reads   ──► RPM ──► bleed filter ──► [neg-control filter] ──► Krona (raw + filtered)
  ├─► [DIAMOND blastx on reads] ──► same filter chain ──► Krona (raw + filtered)
  ├─► [MEGAHIT de novo assembly]
  │     ├─► [racon / medaka polish — Nanopore only]
  │     ├─► [Kraken2 on contigs]    ──► same filter chain ──► Krona pair
  │     └─► [DIAMOND on contigs]    ──► idxstats support filter ──► same filter chain ──► Krona pair
  ├─► [reference-assembly checkpoint] ──► per-(family, accession) consensus FASTA
  └─► MultiQC report (Illumina only)
```

The output of every classifier — reads or contigs — flows through the same cross-sample summary stack: a raw count table, an RPM-normalised table, a bleed-filtered table, and (if you declare negative controls) a Poisson-filtered table. Each classifier also gets a *raw* and a *filtered* Krona HTML per sample so you can compare the unfiltered hits against what survives the cross-sample filters.

## Decision matrix — what to turn on

| If you want…                                                                       | Add…                                                                                                          |
|------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| "What's in my sample?" — fast first look                                           | nothing extra (Kraken2 on reads is the default)                                                                |
| Protein-level sensitivity for divergent viruses                                    | `--run-diamond-reads --diamond-database … --taxids …`                                                          |
| Recover near-full genomes via de novo assembly                                     | `--run-denovo-assembly` (auto-runs Kraken2 on contigs; opt out with `--no-kraken2-contigs`)                    |
| Find viruses missed by Kraken2 but visible via protein homology on assembled contigs | also `--run-diamond-contigs`                                                                                  |
| Automatic per-virus consensus genome from any hit family                            | `--run-reference-assembly --method … --source … --families … --viral-genomes … --viral-taxids …`              |
| Remove host reads first (recommended)                                              | `--deacon-index databases/deacon_indexes/panhuman-1.idx` (or `--host-reference host.fasta`)                    |
| Suppress lab/blank background                                                       | `--negative-controls blank1,blank2`                                                                            |

## Worked example — Illumina, built up incrementally

We will use the same `samples_illumina.csv` from [Setup](setup.md) — the two SARS-CoV-2 amplicon samples `4117` and `61382`. Each step below adds one slice of the pipeline. You can stop at any step that gives you what you need.

### Step A — minimal: Kraken2 on reads

```bash
viralunity meta illumina \
    --sample-sheet     samples_illumina.csv \
    --config-file      results/meta_illumina/config.yml \
    --output           results/meta_illumina/ \
    --run-name         sarscov2 \
    --kraken2-database databases/kraken2 \
    --krona-database   databases/krona/taxonomy \
    --taxdump          databases/taxdump \
    --threads 2 --threads-total 4
```

This is the smallest useful invocation: trim with fastp, run Kraken2 on the trimmed reads, build a Krona plot per sample, and write the cross-sample summary tables. No dehosting, no assembly, no DIAMOND, no reference assembly.

Look at the outputs under `results/meta_illumina/sarscov2/`:

```text
qc/reports/multiqc_report.html
samples/<sample>/
├── fastp.html
├── kraken2_reads.report.txt           # Kraken2 native report
├── kraken2_reads.krona.html           # interactive sunburst — all hits
└── kraken2_reads.filtered.krona.html  # same, but pruned to taxa that pass the bleed filter
metagenomics/taxonomic_assignments/kraken2_reads/
├── results/<sample>.report.txt
├── kraken2_reads_taxa_summary.tsv         # per-rank counts (family/genus/species)
├── kraken2_reads_taxa_summary_RPM.tsv     # adds total_reads + RPM columns
└── kraken2_reads_taxa_summary_RPM.bleed.tsv  # adds bleed_threshold + bleed_pass
```

How to read them:

- **`samples/<sample>/kraken2_reads.report.txt`** — Kraken2's own per-taxon report (5 columns: percent, count, count-direct, rank, taxid, name). Use it as a familiar starting point if you know Kraken2 already.
- **`samples/<sample>/kraken2_reads.krona.html`** — open it in a browser. The interactive sunburst lets you drill from "Viruses" down to specific species. The `.filtered.krona.html` next to it is the same picture pruned to taxa that survived the bleed filter (see [Understanding the filters](#understanding-the-filters) below). Use the raw plot to explore, the filtered plot to act.
- **The three TSVs** progress: raw counts → +RPM → +bleed filter. Each row is one `(sample, tool, mode, rank, taxid)` tuple. Open the `.bleed.tsv` in your spreadsheet of choice and filter to `bleed_pass == True` to get the call list.

```{tip}
For a reference output, compare against `my_results/test_meta_illumina/`. Its config (`my_results/config_meta_illumina.yml`) corresponds to Step D below — the full pipeline.
```

### Step B — add host depletion

Even for non-human samples it is worth dehosting: human DNA from lab personnel and kit components is a common low-level contaminant. We use the pre-built Deacon human index because it is much faster than minimap2 against the full human genome:

```bash
viralunity meta illumina \
    --sample-sheet     samples_illumina.csv \
    --config-file      results/meta_illumina/config.yml \
    --output           results/meta_illumina/ \
    --kraken2-database databases/kraken2 \
    --krona-database   databases/krona/taxonomy \
    --taxdump          databases/taxdump \
    --deacon-index     databases/deacon_indexes/panhuman-1.idx \
    --threads 2 --threads-total 4
```

What changed:

- `--deacon-index` — host depletion runs *before* any classification. New per-sample outputs at `host_filtered/<sample>.R{1,2}.filtered.fastq.gz` and a merged `host_filtered/<sample>.merged.fastq.gz` (the merged file is the denominator used to compute RPM).
- Pass `--host-reference host.fasta` instead if you want minimap2-based depletion against a non-human host you already have as a FASTA. They are mutually exclusive.

### Step C — add de novo assembly + contig classifiers + DIAMOND

```bash
viralunity meta illumina \
    --sample-sheet     samples_illumina.csv \
    --config-file      results/meta_illumina/config.yml \
    --output           results/meta_illumina/ \
    --kraken2-database databases/kraken2 \
    --krona-database   databases/krona/taxonomy \
    --taxdump          databases/taxdump \
    --deacon-index     databases/deacon_indexes/panhuman-1.idx \
    --run-denovo-assembly \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database databases/diamond/viral.dmnd \
    --taxids           databases/diamond/protein2taxid.tsv \
    --threads 4 --threads-total 8
```

What this adds:

- **MEGAHIT** runs on host-filtered reads, producing `denovo_assembly/megahit/<sample>/final.contigs.fa`.
- **Kraken2 on contigs** is on by default whenever `--run-denovo-assembly` is set; the contig-level summary tables land in `metagenomics/taxonomic_assignments/kraken2_contigs/`. (Opt out with `--no-kraken2-contigs`.)
- **DIAMOND blastx** runs on both reads (`diamond_reads/`) and contigs (`diamond_contigs/`). DIAMOND complements Kraken2 because it works at the protein level — useful for divergent viruses whose nucleotide sequence has drifted away from anything in the Kraken2 index.
- The contig-level DIAMOND output is post-filtered to *supported* contigs: reads are mapped back to assembled contigs (minimap2 + samtools idxstats), and a contig is only kept if at least one read maps to it. The output is named `*.diamond.supported.tsv` for this reason.

You now have **four classification streams** running in parallel — `kraken2_reads`, `kraken2_contigs`, `diamond_reads`, `diamond_contigs` — each with its own raw → RPM → bleed (→ neg) summary stack and its own raw / filtered Krona pair under `samples/<sample>/`.

If you want the protein-level signal but only on contigs (and not on reads), drop `--run-diamond-reads`. Conversely, `--no-run-kraken2-reads` turns off the read-level Kraken2 stream.

### Step D — turn on reference assembly

This is the headline feature that makes the meta pipeline more than a classifier. When `--run-reference-assembly` is on, a Snakemake checkpoint reads the filtered taxa tables, picks one or more reference genomes per sample, and runs reference-guided consensus assembly automatically — the same alignment/consensus rules used by `viralunity consensus`.

```bash
viralunity meta illumina \
    --sample-sheet     samples_illumina.csv \
    --config-file      results/meta_illumina/config.yml \
    --output           results/meta_illumina/ \
    --kraken2-database databases/kraken2 \
    --krona-database   databases/krona/taxonomy \
    --taxdump          databases/taxdump \
    --deacon-index     databases/deacon_indexes/panhuman-1.idx \
    --run-denovo-assembly \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database databases/diamond/viral.dmnd \
    --taxids           databases/diamond/protein2taxid.tsv \
    --run-reference-assembly \
    --method  kraken2 \
    --source  reads \
    --families Coronaviridae \
    --reads-count 100 \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids  databases/virus_genomes/genome2taxid.tsv \
    --threads 4 --threads-total 8
```

The new flags:

- `--run-reference-assembly` — turn the feature on.
- `--method` — which classifier's output drives the selection (`kraken2`, `diamond`, or `both`). Required.
- `--source` — which level (`reads`, `contigs`, or `both`). Required.
- `--families` — comma-separated list of viral families to target. Hits to any other family are ignored. Defaults to a broad list of common pathogens; here we restrict to coronaviruses.
- `--reads-count` (or `--contigs-count` for contig sources) — minimum hits to trigger assembly for that family. Below the threshold, no reference is picked.
- `--viral-genomes` + `--viral-taxids` — the genome FASTA and `genome2taxid.tsv` produced by `viralunity get-databases virus-genome`. The pipeline pulls reference sequences from these.

Two new top-level outputs appear in `results/meta_illumina/sarscov2/`:

- **`reference_targets.tsv`** — one row per `(sample, ref_key, accession)` selected by the checkpoint. The `ref_key` is `{family}_{accession}` (e.g. `Coronaviridae_NC_045512`), so two coronaviruses detected in one sample produce two separate ref_keys and two consensus runs.
- **`assembly/<ref_key>/.../consensus/final_consensus/<sample>.consensus.fasta`** — the per-sample consensus against the auto-selected reference. The full output tree under `assembly/<ref_key>/` mirrors what `viralunity consensus` produces (references, BAMs, VCFs, coverage stats).

See [Reference assembly under the hood](#reference-assembly-under-the-hood) below for what the checkpoint actually does.

## Nanopore variant

For Nanopore data the same incremental story applies. The differences from Illumina:

- **No fastp** — long reads skip QC; the variant caller is expected to absorb noise.
- **Polishing is opt-in** — `--run-polish-racon` and/or `--run-polish-medaka --medaka-model r941_min_high_g360` polish the MEGAHIT contigs before contig-level classification. Pick a Medaka model that matches your basecaller.
- **Variant caller in the reference-assembly path** is Clair3. Use `--clair3-model` to pick a model (default `r1041_e82_400bps_sup_v500`).

A full Nanopore run with everything turned on:

```bash
viralunity meta nanopore \
    --sample-sheet     samples_nanopore.csv \
    --config-file      results/meta_nanopore/config.yml \
    --output           results/meta_nanopore/ \
    --kraken2-database databases/kraken2 \
    --krona-database   databases/krona/taxonomy \
    --taxdump          databases/taxdump \
    --deacon-index     databases/deacon_indexes/panhuman-1.idx \
    --run-denovo-assembly \
    --run-polish-medaka --medaka-model r941_min_high_g360 \
    --run-diamond-reads --run-diamond-contigs \
    --diamond-database databases/diamond/viral.dmnd \
    --taxids           databases/diamond/protein2taxid.tsv \
    --run-reference-assembly \
    --method kraken2 --source reads \
    --families Coronaviridae \
    --viral-genomes databases/virus_genomes/viral.genomes.fasta \
    --viral-taxids  databases/virus_genomes/genome2taxid.tsv \
    --clair3-model  r1041_e82_400bps_sup_v500 \
    --threads 4 --threads-total 4
```

Outputs follow the same layout as the Illumina full run; the differences are the missing `qc/reports/multiqc_report.html` and the presence of `denovo_assembly/megahit/<sample>/polished.fasta` (when Medaka polishing is on).

## Understanding the filters

Real metagenomic data is noisy. Without cross-sample filtering, every sample looks like it contains thirty different viruses — most of them index-hopping, kit contamination, or extremely low-level cross-sample bleed during sequencing. ViralUnity's filtering chain is designed to suppress those signals while preserving real positives.

Each `*_taxa_summary*.tsv` file is one step in the chain:

| File                                                    | What it adds                                                                                                  |
|---------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| `<classifier>_taxa_summary.tsv`                         | Raw counts, one row per `(sample, rank, taxid)`.                                                              |
| `<classifier>_taxa_summary_RPM.tsv`                     | Adds `total_reads` (host-filtered read count for the sample) and `rpm` = `count / total_reads × 1e6`.         |
| `<classifier>_taxa_summary_RPM.bleed.tsv`               | Adds `max_rpm`, `bleed_threshold`, `bleed_applied`, `bleed_pass`. (Always produced.)                          |
| `<classifier>_taxa_summary_RPM.bleed.neg.tsv`           | Adds Poisson statistics against negative controls (`mu_bg`, `p_bg`, `neg_pass`). Only when `--negative-controls` was set. |

### RPM normalisation

Different samples have different total read counts. To compare a taxon's signal across samples we normalise by sequencing depth: `rpm = count / total_reads × 1,000,000`, where `total_reads` is the number of reads in the merged host-filtered FASTQ for that sample. The pipeline computes this per sample and joins it onto the summary.

### Bleed filter

For each taxon, the pipeline looks at its maximum RPM across all samples, and sets a threshold at a small fraction of that maximum:

```text
bleed_threshold = max_rpm * bleed_fraction      # bleed_fraction = 0.005 by default
bleed_pass      = rpm >= bleed_threshold
```

For example, if *Coronaviridae* hits 1000 RPM in sample A and 0.2 RPM in sample B:

- `max_rpm` = 1000
- `bleed_threshold` = 1000 × 0.005 = 5 RPM
- Sample A's row passes (`rpm=1000 ≥ 5`).
- Sample B's row fails (`rpm=0.2 < 5`) — flagged as likely cross-sample bleed.

If `max_rpm` is itself very small (`< rpm_floor`, currently 1.0), the filter is a no-op and `bleed_applied` is `False` — there is no reliable signal to filter against, so every row is preserved.

Tune the strictness with `--bleed-fraction` (default `0.005`). Lower values (`0.001`) are stricter; higher values (`0.01`) are more permissive.

### Negative-control filter

If you sequence blanks (no-template controls, extraction blanks, etc.) alongside real samples, declare them as negative controls:

```bash
viralunity meta illumina --negative-controls blank1,blank2 …
```

The sample IDs must match the sample sheet exactly (no `sample-` prefix). For each taxon, the filter computes a background rate from the negative controls:

```text
λ           = total_count_in_negatives / total_reads_in_negatives
μ_bg(sample) = λ × total_reads(sample)
p_bg        = P(Poisson(μ_bg) ≥ observed_count)
neg_pass    = p_bg < negative_p_threshold       # default 0.01
```

In plain English: "given how much of this taxon we see in the blanks, what is the probability of seeing at least as many reads in this sample by chance?" If that probability is below `--negative-p-threshold`, the row passes.

The bleed and negative filters compose: a row appears as a *call* only if it has `bleed_pass == True` **and** (when negative controls were configured) `neg_pass == True`.

### Filtered Krona plots

Every classifier produces two Krona HTMLs per sample — `*.krona.html` (built directly from the raw classifier output) and `*.filtered.krona.html` (pruned to reads/contigs whose lineage has a passing ancestor at family/genus/species). The filtered plot is the visual equivalent of the `*_RPM.bleed[.neg].tsv` call list. Strain-level hits whose species or genus passes are preserved; rows with `taxid==0` are dropped.

## Reference assembly under the hood

When you turn on `--run-reference-assembly`, a Snakemake checkpoint inspects the bleed-filtered TSV for the source you chose (`--method` × `--source` selects one of the four classifier/level combinations) and decides per sample which reference genomes to assemble against. Two strategies are available:

- **`taxid` (default)** — for each passing taxon assigned to a target family, look up the taxid directly in `viral_taxids` (`genome2taxid.tsv`). If the exact taxid is not present (often the classifier picked a strain that is below species in the NCBI taxonomy), resolve to the species-level ancestor via `taxdump` and retry. Fast and deterministic; assumes your classifier database and your `viral_genomes` database came from compatible RefSeq releases.

- **`similarity`** — BLASTN each de novo contig (from MEGAHIT) against `viral_genomes`. For each contig, keep the best-bitscore hit that passes `--blast-qcov` and `--blast-pident`, then validate that its taxid lives under a target family via `taxdump`. Requires `--run-denovo-assembly`, and uses the BLAST index that `viralunity get-databases virus-genome` builds alongside the FASTA. Slower than `taxid` but finds closer relatives when the classifier database and the genome database diverge.

Whichever strategy fires, the selected `(sample, family, accession)` triples are written to `reference_targets.tsv`. From there, each `ref_key = {family}_{accession}` triggers a full reference-guided consensus run, with output keyed under `assembly/<ref_key>/`. When one sample has hits to multiple families, you get multiple `ref_key` directories.

See [Notes — Reference selection strategies](../notes.md#reference-selection-strategies) for the full comparison table and worked invocations.

## Tuning knobs

Sensitivity / specificity of detection:

- **`--bleed-fraction`** (`0.005`) — lower = stricter cross-sample bleed filter.
- **`--negative-p-threshold`** (`0.01`) — lower = stricter blank-background filter.
- **`--minimum-hit-group`** (`4`) — Kraken2's hit-group threshold; raising it makes Kraken2 more conservative.

DIAMOND tuning:

- **`--evalue`** (`1e-10`) — tighter = fewer false homologies.
- **`--diamond-sensitivity`** (`sensitive`) — `mid-sensitive`, `more-sensitive`, `ultra-sensitive` trade compute time for sensitivity to divergent hits.

Reference-assembly triggering:

- **`--reads-count`** (`100`) / **`--contigs-count`** (`1`) — minimum hits to a family before reference assembly fires for that family.
- **`--families`** — restrict reference assembly to a comma-separated list of families.
- **`--reference-selection-strategy`** (`taxid` / `similarity`), **`--blast-qcov`**, **`--blast-pident`** — see [Notes](../notes.md#reference-selection-strategies).

## Common pitfalls

- **Krona conda env must be active.** The `krona` taxonomy must be set up via `viralunity get-databases krona` *while the viralunity conda env is active*, because Krona looks for its taxonomy inside `$CONDA_PREFIX`. If `ktImportTaxonomy` runs with an empty taxonomy you get visually empty HTML plots.
- **`--diamond-database` accepts either `.dmnd` or `.faa`.** If you pass a FASTA, the pipeline runs `diamond makedb` for you on the fly. The output `viral.dmnd` lands next to the FASTA.
- **The `similarity` strategy needs the BLAST index.** `viralunity get-databases virus-genome` builds it alongside `viral.genomes.fasta`. If you brought your own viral genome FASTA, run `makeblastdb -dbtype nucl -in viral.genomes.fasta` manually before running `viralunity meta --reference-selection-strategy similarity`.
- **`--negative-controls` IDs must match the sample sheet exactly** — no `sample-` prefix (that prefix is only used internally in the generated YAML).
- **`--method` and `--source` are required when `--run-reference-assembly` is on.** They do not default; the pipeline will error out if you turn on reference assembly without specifying them.

## Reference

- [Commands reference — `viralunity meta`](../commands.md#viralunity-meta): the complete flag table, including every resource and tuning option.
- [Output layout — Metagenomics pipeline](../output.md#metagenomics-pipeline-viralunity-meta): the full output tree.
- [Notes — Dynamic reference assembly](../notes.md#dynamic-reference-assembly): the full `taxid` vs `similarity` comparison.
- [Notes — Reading the Krona "raw vs filtered" pair](../notes.md): how the filtered Krona is computed from the summary tables.
