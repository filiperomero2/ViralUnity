# ViralUnity — empirical test commands

Two SARS-CoV-2 datasets in `my_test_data/` exercise both the consensus
and metagenomics pipelines (Illumina + Nanopore). Output lands under
`my_results/` which is already covered by the `my_results*` gitignore
rule.

System: 16 cores, 125 GB RAM. `--threads 4 --threads-total 8` lets two
rules run in parallel.

## Setup (once per shell session)

```bash
mamba activate viralunity      # or: conda activate viralunity
cd /home/gevop/projects/ViralUnity

export DATA=/home/gevop/projects/ViralUnity/my_test_data
export RESULTS=/home/gevop/projects/ViralUnity/my_results
export REF=/home/gevop/databases/reference_genomes/nCoV-2019.reference.fasta
export BED=/home/gevop/databases/reference_genomes/nCoV-2019.bed
export ADAPTERS=/home/gevop/databases/reference_genomes/adapters.fa
export DBS=/home/gevop/projects/ViralUnity/databases
export DEACON_INDEX=/home/gevop/databases/deacon_indexes/panhuman-1.k31w15.idx

mkdir -p "$RESULTS"
```

## 1. Sample sheets

Both runs are level-0 (files directly in the data dir, not one subdir
per sample).

```bash
# Illumina: sample = "4117", "61382" (split on "_" before R1/R2)
viralunity create-samplesheet \
    --input "$DATA/illumina_data/" \
    --output "$RESULTS/samples_illumina.csv" \
    --level 0 --separator _ --pattern R1

# Nanopore: sample = "barcode05", "barcode09" (split on "." before .fastq)
viralunity create-samplesheet \
    --input "$DATA/nanopore_data/" \
    --output "$RESULTS/samples_nanopore.csv" \
    --level 0 --separator . --pattern barcode
```

Inspect `samples_illumina.csv` (3-column rows: `sample_id,R1,R2`) and
`samples_nanopore.csv` (2-column rows: `sample_id,fastq`) before
continuing.

## 2. Consensus pipeline

```bash
# --- Illumina ---
viralunity consensus illumina \
    --sample-sheet "$RESULTS/samples_illumina.csv" \
    --config-file "$RESULTS/config_consensus_illumina.yml" \
    --run-name test_consensus_illumina \
    --output "$RESULTS" \
    --reference "$REF" \
    --primer-scheme "$BED" \
    --adapters "$ADAPTERS" \
    --threads 4 --threads-total 8

# --- Nanopore ---
viralunity consensus nanopore \
    --sample-sheet "$RESULTS/samples_nanopore.csv" \
    --config-file "$RESULTS/config_consensus_nanopore.yml" \
    --run-name test_consensus_nanopore \
    --output "$RESULTS" \
    --reference "$REF" \
    --primer-scheme "$BED" \
    --threads 4 --threads-total 8
```

Consensus output appears at
`$RESULTS/test_consensus_illumina/samples/<sample>/consensus.fasta`
(and the Nanopore variant). The cross-sample alignment is at
`assembly/consensus/final_consensus/samples_alignment.fasta`.

## 3. Metagenomics pipeline

Minimal run — Kraken2 on reads only, no dehosting, no de novo
assembly. The Kraken2 viral DB under `databases/kraken2/` is fine for
SARS-CoV-2 detection.

```bash
# --- Illumina ---
viralunity meta illumina \
    --sample-sheet "$RESULTS/samples_illumina.csv" \
    --config-file "$RESULTS/config_meta_illumina.yml" \
    --run-name test_meta_illumina \
    --output "$RESULTS" \
    --kraken2-database "$DBS/kraken2" \
    --krona-database "$DBS/krona/taxonomy" \
    --taxdump "$DBS/taxdump" \
    --adapters "$ADAPTERS" \
    --threads 4 --threads-total 8

# --- Nanopore (no fastp step; no --adapters) ---
viralunity meta nanopore \
    --sample-sheet "$RESULTS/samples_nanopore.csv" \
    --config-file "$RESULTS/config_meta_nanopore.yml" \
    --run-name test_meta_nanopore \
    --output "$RESULTS" \
    --kraken2-database "$DBS/kraken2" \
    --krona-database "$DBS/krona/taxonomy" \
    --taxdump "$DBS/taxdump" \
    --threads 4 --threads-total 8
```

Meta outputs land at
`$RESULTS/test_meta_<datatype>/metagenomics/taxonomic_assignments/kraken2_reads/...`
and per-sample symlinks under `samples/<sample>/`.

## Optional flags worth exercising

Append any of these to the `meta` invocations to test more of the
pipeline surface.

| Goal | Flags |
|------|-------|
| Add DIAMOND on reads | `--run-diamond-reads --diamond-database "$DBS/diamond/viral.dmnd" --taxids "$DBS/diamond/protein2taxid.tsv"` |
| Dehost with Deacon | `--deacon-index "$DEACON_INDEX"` |
| Dehost with minimap2 instead | `--host-reference /home/gevop/databases/reference_genomes/GCF_000001405.40_GRCh38.p14_genomic.fna.gz` |
| De novo assembly + contig classification | `--run-denovo-assembly` (Kraken2 on contigs auto-enabled; add `--run-diamond-contigs` to also run DIAMOND) |
| Dry-run / config-peek (no actual run) | `--create-config-only` |
| Dynamic reference assembly from hits | `--run-reference-assembly --method kraken2 --source reads --families Coronaviridae --reads-count 100 --viral-genomes "$DBS/virus_genomes/viral.genomes.fasta" --viral-taxids "$DBS/virus_genomes/genome2taxid.tsv"` (requires the viral-genomes BLAST index too) |

## Notes

- Pass `--create-config-only` first if you want to inspect the
  generated YAML and see how the new Phase-3 keys
  (`minimap2_consensus_align_flags`, `diamond_max_target_seqs`,
  `kraken2_extra_flags`) appear in the `# parameters` section.
- Disk usage: `/` is at 93% (206 GB free). The consensus runs are
  small; full meta-with-assembly runs on real-world data can write
  tens of GB per sample. For these 4 samples it should be fine, but
  worth watching.
- If a run fails partway through and you re-invoke the same command,
  Snakemake will resume from the last successful step (the meta
  pipeline no longer passes `forceall=True` after the Phase-3
  cleanup).

---

# Segmented viruses

These exercise the `--segmented-reference SEG=PATH` codepath that
ships separate per-segment consensus assemblies. References live
under `/home/gevop/databases/reference_genomes/`. Sample data lives
under `/home/gevop/projects/REVISA/data/`.

Because the REVISA run directories contain many other samples, we
write the CSV sample sheet directly instead of running
`create-samplesheet` with level-1 discovery.

## Guaroa virus — *Orthobunyavirus* (3 segments: S, M, L)

One sample (sample-467), Illumina paired-end shotgun. The reference
files use the `_recent_ref` variant of the segments.

### Setup

```bash
mamba activate viralunity      # or: conda activate viralunity
cd /home/gevop/projects/ViralUnity

export RESULTS=/home/gevop/projects/ViralUnity/my_results
export DBS=/home/gevop/projects/ViralUnity/databases
export REFGENOMES=/home/gevop/databases/reference_genomes
export GUAROA_DATA=/home/gevop/projects/REVISA/data/REVISA32/467_ds.2938356697df4799bbef6c697ac68a7e

mkdir -p "$RESULTS"
```

### 1. Sample sheet (one row, written directly)

CSV column 1 is the bare sample ID; the pipeline prepends `sample-`
when it writes the YAML, so don't pre-prefix here (otherwise you end
up with `sample-sample-467` in every output filename).

```bash
cat > "$RESULTS/samples_guaroa.csv" <<EOF
467,$GUAROA_DATA/467_S10_L001_R1_001.fastq.gz,$GUAROA_DATA/467_S10_L001_R2_001.fastq.gz
EOF
```

### 2. Consensus pipeline (per-segment)

```bash
viralunity consensus illumina \
    --sample-sheet "$RESULTS/samples_guaroa.csv" \
    --config-file "$RESULTS/config_consensus_guaroa.yml" \
    --run-name test_consensus_guaroa \
    --output "$RESULTS" \
    --segmented-reference S=$REFGENOMES/guaroa_virus_segment_S_recent_ref.fasta \
    --segmented-reference M=$REFGENOMES/guaroa_virus_segment_M_recent_ref.fasta \
    --segmented-reference L=$REFGENOMES/guaroa_virus_segment_L_recent_ref.fasta \
    --minimum-coverage 10 \
    --threads 4 --threads-total 8
```

Per-segment consensus FASTAs appear at
`$RESULTS/test_consensus_guaroa/assembly/{S,M,L}/consensus/final_consensus/sample-467.consensus.fasta`,
and a per-segment cross-sample alignment at
`$RESULTS/test_consensus_guaroa/assembly/{S,M,L}/consensus/final_consensus/samples_alignment.fasta`.

### 3. Metagenomics pipeline (classification only)

The meta pipeline doesn't currently have per-segment reference-
assembly logic, so this run is just classification + filters.

```bash
viralunity meta illumina \
    --sample-sheet "$RESULTS/samples_guaroa.csv" \
    --config-file "$RESULTS/config_meta_guaroa.yml" \
    --run-name test_meta_guaroa \
    --output "$RESULTS" \
    --kraken2-database "$DBS/kraken2" \
    --krona-database "$DBS/krona/taxonomy" \
    --taxdump "$DBS/taxdump" \
    --no-kraken2-contigs \
    --threads 4 --threads-total 8
```

`--no-kraken2-contigs` skips contig classification (which would
otherwise require `--run-denovo-assembly`); the default flag combo
fails validation. For a classification-only meta run on targeted
data this is what you want.

## Influenza A H3N2 — *Orthomyxoviridae* (8 segments: S1–S8)

Two samples (sample-1177, sample-1208), Illumina paired-end shotgun.

### Setup (in addition to the env vars above)

```bash
export FLU_DATA=/home/gevop/projects/REVISA/data/REVISA43
```

### 1. Sample sheet (two rows, written directly)

```bash
cat > "$RESULTS/samples_influenza.csv" <<EOF
1177,$FLU_DATA/1177_ds.d9c0f0ef460d45bbb3a7fea896696cf7/1177_S6_L001_R1_001.fastq.gz,$FLU_DATA/1177_ds.d9c0f0ef460d45bbb3a7fea896696cf7/1177_S6_L001_R2_001.fastq.gz
1208,$FLU_DATA/1208_ds.cb05f4155fc64070aa9896a14bc68fcc/1208_S12_L001_R1_001.fastq.gz,$FLU_DATA/1208_ds.cb05f4155fc64070aa9896a14bc68fcc/1208_S12_L001_R2_001.fastq.gz
EOF
```

### 2. Consensus pipeline (per-segment)

```bash
viralunity consensus illumina \
    --sample-sheet "$RESULTS/samples_influenza.csv" \
    --config-file "$RESULTS/config_consensus_influenza.yml" \
    --run-name test_consensus_influenza \
    --output "$RESULTS" \
    --segmented-reference S1=$REFGENOMES/influenza_A_H3N2_segment_1.fasta \
    --segmented-reference S2=$REFGENOMES/influenza_A_H3N2_segment_2.fasta \
    --segmented-reference S3=$REFGENOMES/influenza_A_H3N2_segment_3.fasta \
    --segmented-reference S4=$REFGENOMES/influenza_A_H3N2_segment_4.fasta \
    --segmented-reference S5=$REFGENOMES/influenza_A_H3N2_segment_5.fasta \
    --segmented-reference S6=$REFGENOMES/influenza_A_H3N2_segment_6.fasta \
    --segmented-reference S7=$REFGENOMES/influenza_A_H3N2_segment_7.fasta \
    --segmented-reference S8=$REFGENOMES/influenza_A_H3N2_segment_8.fasta \
    --minimum-coverage 20 \
    --threads 4 --threads-total 8
```

### 3. Metagenomics pipeline (classification only)

```bash
viralunity meta illumina \
    --sample-sheet "$RESULTS/samples_influenza.csv" \
    --config-file "$RESULTS/config_meta_influenza.yml" \
    --run-name test_meta_influenza \
    --output "$RESULTS" \
    --kraken2-database "$DBS/kraken2" \
    --krona-database "$DBS/krona/taxonomy" \
    --taxdump "$DBS/taxdump" \
    --no-kraken2-contigs \
    --threads 4 --threads-total 8
```

## Segmented-virus notes

- The `--segmented-reference` flag is repeatable. Segment labels are
  the keys (`S`, `M`, `L`, `S1`, …); the path is the per-segment
  reference FASTA.
- Outputs are organised under
  `assembly/<segment>/consensus/final_consensus/sample-<id>.consensus.fasta`
  (per-segment, per-sample consensus) and
  `assembly/<segment>/consensus/final_consensus/samples_alignment.fasta`
  (per-segment cross-sample alignment).
- The metagenomics dynamic-reference-assembly path
  (`--run-reference-assembly`) currently extracts a single accession
  per family, so it isn't a faithful test for segmented viruses
  (Influenza/Bunyavirales have one accession per segment). Stick to
  the `consensus` pipeline for per-segment assemblies; use `meta`
  here just for classification/QC.
