# Consensus pipeline tutorial

The consensus pipeline (`viralunity consensus`) turns raw viral reads into a polished consensus genome per sample, aligned to a reference you supply. By the end of this page you will have run it end-to-end on Illumina and Nanopore SARS-CoV-2 data and you will know how to read each output file.

```{note}
This page assumes you have completed [Setup](setup.md): ViralUnity is installed, the `viralunity` conda environment is active, and you have generated `samples_illumina.csv` and `samples_nanopore.csv` from `my_test_data/`.
```

## When to use it

Use `viralunity consensus` when you already know what virus is in the sample and you have a reference genome for it. The pipeline maps reads, calls variants, and produces one consensus FASTA per sample plus a multi-sample alignment ready for tree-building. If you do not know what is in the sample, start with [`viralunity meta`](metagenomics.md) instead.

## What it does

The two data types share the same overall idea — *map → trim primers → call → consensus → align all samples to the reference* — but the variant caller and quality-control steps differ:

```text
Illumina:
  raw FASTQ ──► fastp ──► minimap2 -x sr ──► samtools ampliconclip (BED primers)
            ──► samtools consensus ──► per-base coverage ──► multi-sample alignment
            └► [optional] LoFreq (iSNVs)

Nanopore:
  raw FASTQ ──► (sanitize reference headers) ──► minimap2 -x map-ont
            ──► samtools ampliconclip ──► Clair3 + bcftools (variants + consensus)
            ──► per-base coverage ──► multi-sample alignment
```

You do not pick these tools — they are wired in. You pick:

- the **reference** (`--reference`, or `--segmented-reference SEG=PATH` repeated per segment for influenza-style genomes),
- whether the data is **amplicon** (`--primer-scheme primers.bed` clips amplicon primers) or shotgun,
- the **coverage and allele-frequency thresholds** that decide which bases become `N` (`--minimum-coverage`, `--af-threshold`).

## Worked example 1 — Illumina (SARS-CoV-2 amplicon)

We will assemble the two paired-end samples in `my_test_data/illumina_data/`. The reference and primer scheme below assume the [ARTIC nCoV-2019](https://github.com/artic-network/primer-schemes) bundle, which is the source of `nCoV-2019.reference.fasta` (NCBI `MN908947.3`) and the matching `nCoV-2019.bed`. Substitute paths for whichever copies you have.

### Sample sheet

Already generated in [Setup](setup.md):

```text
4117,my_test_data/illumina_data/4117_S80_L001_R1_001.fastq.gz,my_test_data/illumina_data/4117_S80_L001_R2_001.fastq.gz
61382,my_test_data/illumina_data/61382_S1_L001_R1_001.fastq.gz,my_test_data/illumina_data/61382_S1_L001_R2_001.fastq.gz
```

### Run

```bash
viralunity consensus illumina \
    --sample-sheet  samples_illumina.csv \
    --config-file   results/consensus_illumina/config.yml \
    --output        results/consensus_illumina/ \
    --run-name      sarscov2 \
    --reference     databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --adapters      databases/refs/adapters.fa \
    --minimum-coverage 20 \
    --threads 2 \
    --threads-total 4
```

What each flag does:

- `--sample-sheet` — the CSV from above.
- `--config-file` — where to write the generated Snakemake YAML config (handy to diff against `my_results/config_consensus_illumina.yml`).
- `--output` — base directory for all results.
- `--run-name` — appended to `--output`; everything lands under `results/consensus_illumina/sarscov2/`.
- `--reference` — single-segment reference FASTA. For segmented genomes, use `--segmented-reference` instead (see below).
- `--primer-scheme` — BED file of primer coordinates. Omit it for shotgun data; the primer-clipping step then becomes a pass-through.
- `--adapters` — adapter FASTA passed to fastp. Omit to let fastp auto-detect.
- `--minimum-coverage 20` — positions with fewer than 20 reads after primer clipping become `N` in the consensus.
- `--threads 2 --threads-total 4` — 2 threads per task, 4 cores total across the workflow.

```{tip}
Add `--create-config-only` to write the YAML config and stop. You can then inspect or edit it before running `snakemake --configfile <config> -s viralunity/scripts/consensus_illumina.smk --use-conda -j 4` yourself.
```

Snakemake will print a job table and stream tool output for each rule (`perform_qc`, `map_reads`, `trim_primer_sequences`, `infer_consensus_sequence`, `calculate_assembly_statistics`, `generate_multiqc_report`, `align_consensus_to_reference_genome`, …). On a typical laptop the bundled SARS-CoV-2 samples finish in a few minutes; consult `results/consensus_illumina/sarscov2/benchmark.tsv` after the run for per-rule timing.

### Tour the outputs

After the run, the relevant paths under `results/consensus_illumina/sarscov2/` are:

```text
samples/<sample>/
├── consensus.fasta            # final consensus sequence
├── consensus.vcf.gz           # variants relative to the reference
├── fastp.html                 # fastp QC report
├── raw_mapped_reads.bam       # post-mapping, before primer clipping
├── trimmed_mapped_reads.bam   # primer-clipped — the BAM used for consensus
└── stats_summary.csv          # per-sample mapped/coverage stats
assembly/
├── assembly_stats_summary.csv          # the per-sample stats combined into one CSV
├── coverage_stats/<sample>.table_cov_basewise.txt  # per-base depth
└── consensus/final_consensus/
    └── samples_alignment.fasta          # all samples + reference, MSA-ready
qc/reports/multiqc_report.html           # combined fastp/QC report
benchmark.tsv                            # wall time + memory per rule per sample
```

How to read each one:

**`samples/<sample>/consensus.fasta`** — your finished genome. Long runs of `N` indicate stretches with coverage below `--minimum-coverage` (or where every read disagreed with the reference but no allele exceeded `--af-threshold`). A first sanity check is the proportion of non-N bases — `assembly_stats_summary.csv` reports it as `horizontal_coverage`.

**`samples/<sample>/consensus.vcf.gz`** — the differences between your sample and the reference, called from the consensus FASTA via GSAlign. Inspect with:

```bash
bcftools view results/consensus_illumina/sarscov2/samples/4117/consensus.vcf.gz | head -30
```

**`assembly/coverage_stats/<sample>.table_cov_basewise.txt`** — three columns: `RNAME`, `POS`, `DEPTH`. Useful for spotting drop-outs:

```bash
awk '$3 < 20' results/consensus_illumina/sarscov2/assembly/coverage_stats/4117.table_cov_basewise.txt | head
```

**`samples/<sample>/{raw,trimmed}_mapped_reads.bam`** — both exist deliberately. `raw_mapped_reads.bam` is what minimap2 produced; `trimmed_mapped_reads.bam` is the same BAM after `samtools ampliconclip` removed primer sequences (only different when you passed `--primer-scheme`). The consensus is called from the trimmed BAM.

**`assembly/assembly_stats_summary.csv`** — one row per sample, with `number_of_reads`, `number_of_trim_paired_reads`, `number_of_mapped_reads`, `average_depth`, `percentage_above_{10,100,1000}x`, and `horizontal_coverage`. Quick way to spot low-coverage or poorly-mapping samples without opening each BAM.

**`assembly/consensus/final_consensus/samples_alignment.fasta`** — all per-sample consensuses plus the reference, aligned (built by `minimap2` followed by `gofasta sam toMultiAlign`). Drop this straight into a tree-builder such as IQ-TREE for a quick phylogeny.

**`benchmark.tsv`** — every Snakemake task's runtime, memory, and CPU. Useful when you scale up to a real run.

```{tip}
For a reference run, compare your output tree against `my_results/test_consensus_illumina/` in this checkout — it was produced by running exactly this command on the same input.
```

### Optional: intra-host SNVs

If you care about within-sample variation (e.g., looking for sub-consensus mutations), add `--run-isnv`:

```bash
viralunity consensus illumina \
    --sample-sheet samples_illumina.csv \
    --config-file  results/consensus_illumina_isnv/config.yml \
    --output       results/consensus_illumina_isnv/ \
    --reference    databases/refs/nCoV-2019.reference.fasta \
    --primer-scheme databases/refs/nCoV-2019.bed \
    --run-isnv \
    --af-isnv-threshold 0.05 \
    --threads 2 --threads-total 4
```

This adds a LoFreq pass on the primer-clipped BAM. The output VCFs land under `assembly/isnvs/<sample>.isnvs.vcf.gz`, and a per-sample tally is written to `isnvs/isnvs_summary.tsv`. `--af-isnv-threshold` sets the minimum allele frequency reported (LoFreq filters consensus-level variants out so this captures the sub-consensus range).

## Worked example 2 — Nanopore (SARS-CoV-2)

Same virus, different chemistry. The reads are in `my_test_data/nanopore_data/`. Nanopore's higher per-base error rate calls for a deep-learning variant caller (Clair3) and a stricter default allele-frequency threshold than Illumina.

### Sample sheet (2 columns)

```text
barcode05,my_test_data/nanopore_data/barcode05.fastq
barcode09,my_test_data/nanopore_data/barcode09.fastq
```

### Run

```bash
viralunity consensus nanopore \
    --sample-sheet samples_nanopore.csv \
    --config-file  results/consensus_nanopore/config.yml \
    --output       results/consensus_nanopore/ \
    --run-name     sarscov2 \
    --reference    databases/refs/nCoV-2019.reference.fasta \
    --clair3-model r1041_e82_400bps_sup_v500 \
    --minimum-coverage 20 \
    --minimum-map-quality 30 \
    --threads 4 --threads-total 4
```

The interesting differences from the Illumina invocation:

- No `--adapters` and no `fastp` step — the workflow skips QC entirely for Nanopore. The variant caller is expected to absorb noisy bases.
- `--clair3-model` picks the appropriate Clair3 model for your basecaller. The default `r1041_e82_400bps_sup_v500` matches recent (R10.4.1) SUP-basecalled reads; pick a different one if your data was basecalled differently. The list of model names is in the [Clair3 model zoo](https://github.com/HKU-BAL/Clair3#pre-trained-models).
- `--minimum-map-quality 30` filters reads with MAPQ below 30 before variant calling. Tighten this for very noisy runs.
- The pipeline silently sanitizes the reference's FASTA headers (replacing `/`, `|`, `,`, `~`, and spaces with `_`) before use — clair3 turns the seq ID into a directory name, so the sanitization avoids cryptic filesystem errors. The sanitized copy lands at `results/consensus_nanopore/sarscov2/reference/reference.sanitized.fasta`.

### What is different in the output

Most of the output layout matches the Illumina run, with a couple of additions and substitutions:

- `reference/reference.sanitized.fasta` — the sanitized reference used by every downstream rule. If you need the exact contig names that appear in the BAM and VCF, look here, not at your input FASTA.
- `samples/<sample>/consensus.vcf.gz` — produced by Clair3 + `bcftools` rather than from the consensus directly. You will see `QUAL` values reflective of clair3's confidence.
- No `fastp.html` per sample, no `multiqc_report.html`.

The headline `consensus.fasta` and `samples_alignment.fasta` files behave identically to the Illumina case.

```{tip}
The default `--af-threshold` for Nanopore is intentionally lower-leaning than what you might expect from Illumina, because Clair3 already filters by its own quality score. If you see too many heterozygous-looking sites in your consensus, raise it: `--af-threshold 0.7` or higher.
```

## Segmented viruses

For multi-segment genomes (influenza, the bunyaviruses, etc.) repeat `--segmented-reference SEG=PATH` once per segment instead of using `--reference`:

```bash
viralunity consensus illumina \
    --sample-sheet samples_flu.csv \
    --config-file  results/flu/config.yml \
    --output       results/flu/ \
    --segmented-reference HA=refs/flu/HA.fasta \
    --segmented-reference NA=refs/flu/NA.fasta \
    --segmented-reference PB1=refs/flu/PB1.fasta \
    # … one per segment …
    --threads 4 --threads-total 8
```

The workflow dispatches a specialised segmented variant that processes each segment in parallel. Outputs are keyed by segment under `assembly/<segment>/consensus/final_consensus/` and per-sample symlinks land under `samples/<sample>/<segment>/`. A real config example is at `my_results/config_consensus_influenza.yml` (8 segments) and `my_results/config_consensus_guaroa.yml` (3 segments).

See [Notes — Segmented viruses](../notes.md#segmented-viruses) for the full discussion.

## Tuning consensus quality

The three knobs you will reach for most often:

- **`--minimum-coverage`** (default `20`) — any reference position with fewer than this many reads after primer clipping becomes `N`. Lower it (e.g. `10`) for low-depth samples where you would rather see a tentative base than a hole; raise it for high-confidence assemblies.
- **`--af-threshold`** (default `0.51`) — minimum allele frequency to call a variant into the consensus. At `0.51` the consensus is the majority allele. Raise it (e.g. `0.7`) for noisier data; lower it to capture ambiguity codes.
- **`--minimum-read-length`** (default `50`) — reads below this length are dropped at QC time on Illumina, and used as a minimum-length filter on primer-clipped reads on Nanopore.

Nanopore has a few extra knobs:

- **`--variant-quality`** (default `20`) — Clair3 minimum QUAL.
- **`--variant-depth`** (default `10`) — minimum alt-allele read support.
- **`--chunk-size`** (default `10000`) — Clair3 chunk size; tune only if Clair3 runs out of memory on huge references.

The full set is in the [Commands reference](../commands.md#viralunity-consensus).

## Common pitfalls

- **Sample sheet column count is the only check.** A Nanopore CSV that accidentally has three columns is silently parsed as Illumina (and vice versa). Always confirm the row count and column count before running.
- **Reference header sanitization is silent.** On Nanopore, if your downstream tooling expects the original reference contig names, read them from `reference/reference.sanitized.fasta`, not your input FASTA.
- **Primer scheme contig names must match the reference.** The BED file must use the same chromosome/contig names as the reference FASTA. If `samtools ampliconclip` reports zero primers clipped, you have a name mismatch.
- **`--reference` and `--segmented-reference` are mutually exclusive.** Pass one or the other.

## Reference

- [Commands reference — `viralunity consensus`](../commands.md#viralunity-consensus): the complete flag table.
- [Output layout — Consensus pipeline](../output.md#consensus-pipeline-viralunity-consensus): the full output tree.
- [Notes — Segmented viruses](../notes.md#segmented-viruses): more on multi-segment runs and the segmented output layout.

Next: [Metagenomics pipeline tutorial →](metagenomics.md)
