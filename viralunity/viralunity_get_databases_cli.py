"""Click CLI for viralunity get-databases command group."""

import json
import logging
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import click
from Bio import SeqIO

logger = logging.getLogger(__name__)


def _make_db_dir(path: str, db_name: str) -> Path:
    """Create and return the database subdirectory."""
    db_dir = Path(path) / db_name
    db_dir.mkdir(parents=True, exist_ok=True)
    click.echo(f"Database directory: {db_dir}")
    return db_dir


def _run(cmd: list, cwd: str | None = None) -> None:
    """Run a shell command, streaming output and raising on failure."""
    click.echo(f"$ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise click.ClickException(
            f"Command failed with exit code {result.returncode}: {' '.join(str(c) for c in cmd)}"
        )


@click.group(name="get-databases")
def get_databases():
    """Download and set up reference databases for ViralUnity pipelines."""


@get_databases.command("kraken2")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the kraken2/ subdirectory will be created.",
)
@click.option(
    "--url",
    default="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz",
    show_default=True,
    help="URL of the Kraken2 pre-built index to download.",
)
def get_kraken2(path, url):
    """Download and extract a Kraken2 pre-built index (default: k2_viral).

    Creates {path}/kraken2/ and downloads the index archive there.
    Check https://benlangmead.github.io/aws-indexes/k2 for available indexes.
    """
    db_dir = _make_db_dir(path, "kraken2")
    archive_name = url.split("/")[-1]
    archive_path = db_dir / archive_name

    click.echo(f"Downloading Kraken2 database from:\n  {url}")
    _run(["wget", "-c", "-O", str(archive_path), url])

    click.echo("Extracting archive...")
    _run(["tar", "-xzvf", str(archive_path), "-C", str(db_dir)])

    click.echo(f"\nKraken2 database ready at: {db_dir}")
    click.echo(f"Use --kraken2-database {db_dir} in your viralunity meta commands.")


@get_databases.command("krona")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the krona/taxonomy/ subdirectory will be created.",
)
def get_krona(path):
    """Set up the Krona taxonomy database.

    Creates {path}/krona/taxonomy/, removes the default taxonomy bundled with
    the conda environment, symlinks the new directory, and runs
    ktUpdateTaxonomy.sh to populate it.

    Requires the viralunity conda environment to be active (CONDA_PREFIX must
    be set).
    """
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        raise click.ClickException(
            "CONDA_PREFIX is not set. Please activate the viralunity conda "
            "environment before running this command."
        )

    taxonomy_dir = _make_db_dir(path, "krona") / "taxonomy"
    taxonomy_dir.mkdir(parents=True, exist_ok=True)

    conda_krona_taxonomy = Path(conda_prefix) / "opt" / "krona" / "taxonomy"

    if conda_krona_taxonomy.exists() or conda_krona_taxonomy.is_symlink():
        click.echo(f"Removing existing Krona taxonomy: {conda_krona_taxonomy}")
        if conda_krona_taxonomy.is_symlink():
            conda_krona_taxonomy.unlink()
        else:
            shutil.rmtree(conda_krona_taxonomy)

    click.echo(f"Symlinking {conda_krona_taxonomy} - {taxonomy_dir}")
    conda_krona_taxonomy.symlink_to(taxonomy_dir.resolve())

    click.echo("Running ktUpdateTaxonomy.sh ...")
    _run(["ktUpdateTaxonomy.sh", str(taxonomy_dir)])

    click.echo(f"\nKrona taxonomy ready at: {taxonomy_dir}")
    click.echo(f"Use --krona-database {taxonomy_dir} in your viralunity meta commands.")


@get_databases.command("taxdump")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the taxdump/ subdirectory will be created.",
)
@click.option(
    "--url",
    default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    show_default=True,
    help="URL of the NCBI taxdump archive.",
)
def get_taxdump(path, url):
    """Download and extract the NCBI taxdump (nodes.dmp, names.dmp, ...).

    Creates {path}/taxdump/ and extracts the archive there.
    Pass {path}/taxdump/ to --taxdump in viralunity meta commands.
    """
    db_dir = _make_db_dir(path, "taxdump")
    archive_name = url.split("/")[-1]
    archive_path = db_dir / archive_name

    click.echo(f"Downloading NCBI taxdump from:\n  {url}")
    _run(["wget", "-c", "-O", str(archive_path), url])

    click.echo("Extracting archive...")
    _run(["tar", "-xzvf", str(archive_path), "-C", str(db_dir)])

    for required in ("nodes.dmp", "names.dmp"):
        if not (db_dir / required).exists():
            raise click.ClickException(
                f"Expected file not found after extraction: {db_dir / required}"
            )

    click.echo(f"\nTaxdump ready at: {db_dir}")
    click.echo(f"Use --taxdump {db_dir} in your viralunity meta commands.")


# Characters that, on their own, would make a sequence look like DNA/RNA. All
# six are also valid amino-acid one-letter codes (A, C, G, T, U, N), so this
# alphabet is only meaningful in combination with a minimum length check.
_DNA_LIKE_ALPHABET = set("ACGTUN")


def _looks_like_dna(seq: str, min_length: int = 20) -> bool:
    """Return ``True`` if ``seq`` looks like a DNA/RNA nucleotide sequence.

    The check is intentionally conservative: a sequence is considered DNA only
    if it is at least ``min_length`` residues long *and* every non-gap, non-
    whitespace character (after upper-casing) belongs to the DNA/RNA alphabet
    ``{A, C, G, T, U, N}``. Real proteins of length >= 20 almost always contain
    at least one amino-acid letter outside this set (M, L, K, R, ...), so
    false positives are negligible in practice.
    """
    cleaned = "".join(seq.split()).replace("-", "").replace("*", "").upper()
    if len(cleaned) < min_length:
        return False
    return all(c in _DNA_LIKE_ALPHABET for c in cleaned)


def _parse_data_report(report_path: Path) -> dict:
    org2taxid = {}
    with open(report_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            org_data = record.get("virus", record.get("organism", {}))
            org_name = org_data.get("organismName", "")
            taxid = str(org_data.get("taxId", ""))
            if org_name and taxid:
                org2taxid[org_name] = taxid
    return org2taxid


def _reformat_protein_fasta(
    protein_files: list,
    output_faa: Path,
    org2taxid: dict,
) -> dict:
    """Extract protein_acc -> taxid map.

    NCBI Datasets newer format headers:
      >NP_056776.2:1-3391 polyprotein [organism=dengue virus type 2]

    This function extracts the organism name to map to the taxID from data_report.jsonl,
    and returns a mapping of the protein accession (the first token after >) to the taxID.

    Records whose sequence is entirely DNA/RNA are skipped. NCBI Datasets has
    been observed to leak nucleotide CDS sequences into ``protein.faa`` for
    some virus genomes; those records use a genome accession (``NC_*``) as
    their header instead of a protein accession (``NP_*`` / ``YP_*``) and
    would otherwise cause ``diamond makedb`` to abort with::

        Error: The sequences are expected to be proteins but only contain DNA letters.
    """
    organism_re = re.compile(r"\[organism=([^\]]+)\]")
    taxid_map = {}  # protein_acc -> taxid (written to protein2taxid.tsv)
    written = 0
    skipped_dna = 0
    skipped_dna_examples = []

    with open(output_faa, "w") as out_handle:
        for prot_file in protein_files:
            for record in SeqIO.parse(prot_file, "fasta"):
                if _looks_like_dna(str(record.seq)):
                    skipped_dna += 1
                    if len(skipped_dna_examples) < 5:
                        skipped_dna_examples.append(record.id)
                    continue

                org_m = organism_re.search(record.description)
                if org_m:
                    org_name = org_m.group(1)
                    taxid = org2taxid.get(org_name, "0")

                    # Extract protein accession (everything before ':'), removing range
                    protein_acc = record.id.split(":")[0]
                    taxid_map[protein_acc] = taxid
                    record.id = protein_acc
                    record.description = ""
                    written += 1

                SeqIO.write(record, out_handle, "fasta")

    click.echo(f"  Processed {written} protein sequences.")
    if skipped_dna:
        preview = ", ".join(skipped_dna_examples)
        suffix = ", ..." if skipped_dna > len(skipped_dna_examples) else ""
        click.echo(
            f"  WARNING: Skipped {skipped_dna} DNA/RNA record(s) found in "
            f"protein.faa (e.g. {preview}{suffix}). "
            "NCBI Datasets sometimes leaks nucleotide CDS sequences into the "
            "protein output; those records were removed so diamond makedb "
            "can build the database."
        )
    return taxid_map


@get_databases.command("diamond")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the diamond/ subdirectory will be created.",
)
@click.option(
    "--taxon",
    default="Viruses",
    show_default=True,
    help="NCBI taxon name to download (e.g. 'Viruses', 'coronaviridae').",
)
@click.option(
    "--refseq/--no-refseq",
    default=True,
    show_default=True,
    help="Limit to RefSeq genomes only.",
)
@click.option(
    "--threads",
    default=1,
    show_default=True,
    type=int,
    help="Threads for diamond makedb.",
)
@click.option(
    "--skip-makedb",
    is_flag=True,
    default=False,
    help="Download and reformat files only; skip diamond makedb.",
)
def get_diamond(path, taxon, refseq, threads, skip_makedb):
    """Download viral proteins via NCBI Datasets and build a Diamond database.

    Uses 'datasets download virus genome' with --include protein to fetch
    RefSeq viral protein sequences.  The downloaded proteins are reformatted
    so that each FASTA header becomes:

        >genome_accession

    A two-column protein2taxid.tsv mapping file (genome_accession -> TaxID) is
    derived from the bundled data_report.jsonl and saved alongside the database.
    """
    db_dir = _make_db_dir(path, "diamond")
    raw_dir = db_dir / "_ncbi_download"
    raw_dir.mkdir(parents=True, exist_ok=True)
    zip_path = raw_dir / "ncbi_dataset.zip"

    cmd = [
        "datasets",
        "download",
        "virus",
        "genome",
        "taxon",
        taxon,
        "--include",
        "protein",
        "--annotated",
        "--filename",
        str(zip_path),
    ]
    if refseq:
        cmd.append("--refseq")

    click.echo(f"Downloading viral proteins for taxon '{taxon}' from NCBI Datasets...")
    _run(cmd)

    click.echo("Extracting archive...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(raw_dir)

    report_candidates = list(raw_dir.rglob("data_report.jsonl"))
    if not report_candidates:
        raise click.ClickException(
            "data_report.jsonl not found in downloaded package. "
            "The NCBI Datasets package structure may have changed."
        )
    report_path = report_candidates[0]
    click.echo(f"Parsing taxonomy metadata from: {report_path}")
    acc2taxid = _parse_data_report(report_path)
    click.echo(f"  Found {len(acc2taxid)} genome accessions with TaxIDs.")

    protein_files = list(raw_dir.rglob("protein.faa"))
    if not protein_files:
        raise click.ClickException(
            "No protein.faa files found in downloaded package. "
            "Ensure --annotated genomes with proteins are available for this taxon."
        )
    click.echo(f"Found {len(protein_files)} protein file(s).")

    reformatted_faa = db_dir / "viral.protein.faa"
    click.echo(f"Reformatting headers -> {reformatted_faa}")
    taxid_map = _reformat_protein_fasta(protein_files, reformatted_faa, acc2taxid)

    taxid_map_path = db_dir / "protein2taxid.tsv"
    with open(taxid_map_path, "w") as f:
        for genome_acc, taxid in sorted(taxid_map.items()):
            f.write(f"{genome_acc}\t{taxid}\n")
    click.echo(f"Taxonomy mapping written to: {taxid_map_path}")
    click.echo(f"  {len(taxid_map)} protein accessions mapped.")

    shutil.rmtree(raw_dir)
    click.echo("Cleaned up temporary download files.")

    if not skip_makedb:
        dmnd = db_dir / "viral.dmnd"
        click.echo(f"Building Diamond database ({threads} thread(s))...")
        _run(
            [
                "diamond",
                "makedb",
                "--in",
                str(reformatted_faa),
                "--db",
                str(dmnd),
                "--threads",
                str(threads),
            ]
        )
        click.echo("\nDiamond database ready.")
        click.echo(f"  DB:          {dmnd}")
        click.echo(f"  Taxid map:   {taxid_map_path}")
        click.echo("\nUse in viralunity meta commands:")
        click.echo(f"  --diamond-database  {dmnd}")
        click.echo(f"  --taxids  {taxid_map_path}")
    else:
        click.echo(
            f"\nFiles ready in {db_dir} (--skip-makedb was set).\n"
            "Run 'diamond makedb --in viral.protein.faa --db viral.dmnd' manually,\n"
            "then use:\n"
            f"  --diamond-database  {db_dir}/viral.dmnd\n"
            f"  --taxids  {taxid_map_path}"
        )


def _parse_genome_data_report(report_path: Path) -> dict:
    """Parse data_report.jsonl and return accession -> taxid mapping.

    Each JSONL record is expected to have an 'accession' field and a
    'virus'/'organism' block containing 'taxId'.
    """
    acc2taxid = {}
    with open(report_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            accession = record.get("accession", "")
            org_data = record.get("virus", record.get("organism", {}))
            taxid = str(org_data.get("taxId", ""))
            if accession and taxid:
                acc2taxid[accession] = taxid
    return acc2taxid


def _reformat_genome_fasta(
    genome_files: list,
    output_fasta: Path,
    acc2taxid: dict,
) -> dict:
    """Reformat genome FASTA to use accession-only headers.

    Returns a mapping of accession -> taxid for all sequences written.
    """
    taxid_map = {}
    written = 0

    with open(output_fasta, "w") as out_handle:
        for genome_file in genome_files:
            for record in SeqIO.parse(genome_file, "fasta"):
                # The accession is the first token in the header (e.g. NC_045512.2)
                accession = record.id
                taxid = acc2taxid.get(accession, "0")
                taxid_map[accession] = taxid
                record.id = accession
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")
                written += 1

    click.echo(f"  Processed {written} genome sequences.")
    return taxid_map


@get_databases.command("clean-protein-fasta")
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Input protein FASTA file (e.g. databases/diamond/viral.protein.faa).",
)
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(dir_okay=False, writable=True),
    help=(
        "Output FASTA path. Defaults to '<input>.cleaned.faa'. If equal to "
        "--input, the original file is replaced atomically (a backup is kept "
        "at '<input>.with_dna.bak' unless --no-backup is given)."
    ),
)
@click.option(
    "--no-backup",
    is_flag=True,
    default=False,
    help="When replacing the input file in place, do not keep a .with_dna.bak backup.",
)
@click.option(
    "--min-dna-length",
    default=20,
    show_default=True,
    type=int,
    help="Minimum sequence length to classify a record as DNA (smaller records are kept).",
)
def clean_protein_fasta(input_path, output_path, no_backup, min_dna_length):
    """Strip nucleotide records from a protein FASTA file.

    NCBI Datasets occasionally leaks DNA CDS sequences into the ``protein.faa``
    bundled with virus genome packages (the offending records use a genome
    accession such as ``NC_xxxxxx.1`` as their header). ``diamond makedb`` then
    refuses to build the database with::

        Error: The sequences are expected to be proteins but only contain DNA letters.

    This command rewrites the file, dropping any record whose sequence is at
    least ``--min-dna-length`` residues long and made up exclusively of the
    DNA/RNA alphabet ``{A, C, G, T, U, N}``. All other records are passed
    through untouched.
    """
    input_p = Path(input_path)
    if output_path is None:
        output_p = (
            input_p.with_suffix(input_p.suffix + ".cleaned.faa")
            if input_p.suffix
            else Path(str(input_p) + ".cleaned.faa")
        )
    else:
        output_p = Path(output_path)

    in_place = output_p.resolve() == input_p.resolve()
    tmp_p = output_p.with_suffix(output_p.suffix + ".tmp") if in_place else output_p

    kept = 0
    dropped = 0
    dropped_examples = []

    click.echo(f"Cleaning protein FASTA: {input_p}")
    with open(tmp_p, "w") as out_handle:
        for record in SeqIO.parse(str(input_p), "fasta"):
            if _looks_like_dna(str(record.seq), min_length=min_dna_length):
                dropped += 1
                if len(dropped_examples) < 5:
                    dropped_examples.append(record.id)
                continue
            SeqIO.write(record, out_handle, "fasta")
            kept += 1

    if in_place:
        if not no_backup:
            backup_p = input_p.with_suffix(input_p.suffix + ".with_dna.bak")
            shutil.move(str(input_p), str(backup_p))
            click.echo(f"Backed up original to: {backup_p}")
        else:
            input_p.unlink()
        shutil.move(str(tmp_p), str(output_p))

    click.echo(f"Kept {kept} protein record(s).")
    if dropped:
        preview = ", ".join(dropped_examples)
        suffix = ", ..." if dropped > len(dropped_examples) else ""
        click.echo(f"Dropped {dropped} DNA/RNA record(s) (e.g. {preview}{suffix}).")
    else:
        click.echo("No DNA/RNA records found.")
    click.echo(f"Cleaned FASTA written to: {output_p}")


@get_databases.command("virus-genome")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the virus_genomes/ subdirectory will be created.",
)
@click.option(
    "--taxon",
    default="Viruses",
    show_default=True,
    help="NCBI taxon name to download (e.g. 'Viruses', 'coronaviridae').",
)
@click.option(
    "--refseq/--no-refseq",
    default=True,
    show_default=True,
    help="Limit to RefSeq genomes only.",
)
@click.option(
    "--skip-makeblastdb",
    is_flag=True,
    default=False,
    help="Download and reformat files only; skip makeblastdb index creation.",
)
def get_virus_genome(path, taxon, refseq, skip_makeblastdb):
    """Download viral genomes via NCBI Datasets and build a BLAST database.

    Uses 'datasets download virus genome' with --include genome to fetch
    viral genome sequences.  The downloaded genomes are reformatted so
    that each FASTA header becomes:

        >accession

    A two-column genome2taxid.tsv mapping file (accession -> TaxID) is
    derived from the bundled data_report.jsonl and saved alongside the
    FASTA file.  A nucleotide BLAST database is then built with makeblastdb
    for use with the reference assembly similarity strategy.
    """
    db_dir = _make_db_dir(path, "virus_genomes")
    raw_dir = db_dir / "_ncbi_download"
    raw_dir.mkdir(parents=True, exist_ok=True)
    zip_path = raw_dir / "ncbi_dataset.zip"

    cmd = [
        "datasets",
        "download",
        "virus",
        "genome",
        "taxon",
        taxon,
        "--include",
        "genome",
        "--filename",
        str(zip_path),
    ]
    if refseq:
        cmd.append("--refseq")

    click.echo(f"Downloading viral genomes for taxon '{taxon}' from NCBI Datasets...")
    _run(cmd)

    click.echo("Extracting archive...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(raw_dir)

    # Parse taxonomy metadata
    report_candidates = list(raw_dir.rglob("data_report.jsonl"))
    if not report_candidates:
        raise click.ClickException(
            "data_report.jsonl not found in downloaded package. "
            "The NCBI Datasets package structure may have changed."
        )
    report_path = report_candidates[0]
    click.echo(f"Parsing taxonomy metadata from: {report_path}")
    acc2taxid = _parse_genome_data_report(report_path)
    click.echo(f"  Found {len(acc2taxid)} genome accessions with TaxIDs.")

    # Find genome FASTA files
    genome_files = list(raw_dir.rglob("*.fna"))
    if not genome_files:
        raise click.ClickException("No genome FASTA files (.fna) found in downloaded package.")
    click.echo(f"Found {len(genome_files)} genome file(s).")

    # Reformat FASTA
    reformatted_fasta = db_dir / "viral.genomes.fasta"
    click.echo(f"Reformatting headers -> {reformatted_fasta}")
    taxid_map = _reformat_genome_fasta(genome_files, reformatted_fasta, acc2taxid)

    # Write genome2taxid.tsv
    taxid_map_path = db_dir / "genome2taxid.tsv"
    with open(taxid_map_path, "w") as f:
        for acc, taxid in sorted(taxid_map.items()):
            f.write(f"{acc}\t{taxid}\n")
    click.echo(f"Taxonomy mapping written to: {taxid_map_path}")
    click.echo(f"  {len(taxid_map)} genome accessions mapped.")

    # Clean up
    shutil.rmtree(raw_dir)
    click.echo("Cleaned up temporary download files.")

    if not skip_makeblastdb:
        click.echo("Building nucleotide BLAST database (makeblastdb)...")
        _run(
            [
                "makeblastdb",
                "-in",
                str(reformatted_fasta),
                "-dbtype",
                "nucl",
            ]
        )
        click.echo("BLAST database index files created alongside the FASTA.")

    click.echo("\nViral genome database ready.")
    click.echo(f"  FASTA:       {reformatted_fasta}")
    click.echo(f"  Taxid map:   {taxid_map_path}")
    if not skip_makeblastdb:
        click.echo(f"  BLAST DB:    {reformatted_fasta} (index files: .nhr/.nin/.nsq/...)")
    click.echo("\nUse in viralunity meta commands:")
    click.echo(f"  --viral-genomes  {reformatted_fasta}")
    click.echo(f"  --viral-taxids   {taxid_map_path}")


@get_databases.command("host-genome")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the host_genomes/ subdirectory will be created.",
)
@click.option(
    "--accession",
    required=True,
    help="NCBI genome accession ID (e.g. GCA_000001405.29).",
)
def get_host_genome(path, accession):
    """Download a host genome using NCBI Datasets.

    Creates {path}/host_genomes/ and downloads the genome
    sequence and metadata using the provided accession ID.

    Requires 'datasets' (NCBI Datasets CLI) to be installed and in PATH.
    """
    db_dir = Path(path) / "host_genomes"
    db_dir.mkdir(parents=True, exist_ok=True)

    raw_dir = db_dir / f"_ncbi_download_{accession}"
    raw_dir.mkdir(parents=True, exist_ok=True)
    zip_path = raw_dir / "ncbi_dataset.zip"

    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--filename",
        str(zip_path),
    ]

    click.echo(f"Downloading host genome {accession} from NCBI Datasets...")
    _run(cmd)

    click.echo("Extracting archive...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(raw_dir)

    fasta_candidates = list(raw_dir.rglob("*.fna"))
    if not fasta_candidates:
        shutil.rmtree(raw_dir)
        raise click.ClickException(f"FASTA file not found for accession {accession}.")

    fasta_src = fasta_candidates[0]
    fasta_dest = db_dir / f"{accession}.fasta"
    shutil.move(str(fasta_src), str(fasta_dest))

    report_candidates = list(raw_dir.rglob("assembly_data_report.jsonl"))
    info_dest = db_dir / f"{accession}.info.txt"

    if report_candidates:
        report_path = report_candidates[0]
        with open(report_path) as f:
            for line in f:
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                    with open(info_dest, "w") as out:
                        out.write(f"Accession: {accession}\n")
                        org = record.get("organism", {})
                        out.write(f"Organism Name: {org.get('organismName', 'Unknown')}\n")
                        out.write(f"TaxID: {org.get('taxId', 'Unknown')}\n")
                        asm = record.get("assemblyInfo", {})
                        out.write(f"Assembly Name: {asm.get('assemblyName', 'Unknown')}\n")
                        out.write(f"Assembly Level: {asm.get('assemblyLevel', 'Unknown')}\n")
                        binfo = asm.get("biosample", {})
                        if binfo:
                            out.write(f"BioSample: {binfo.get('accession', 'Unknown')}\n")

                        out.write("\n--- Full JSON Report ---\n")
                        out.write(json.dumps(record, indent=2))
                    break
                except json.JSONDecodeError:
                    pass

    shutil.rmtree(raw_dir)
    click.echo("\nHost genome downloaded successfully.")
    click.echo(f"  FASTA: {fasta_dest}")
    if info_dest.exists():
        click.echo(f"  Info:  {info_dest}")


@get_databases.command("deacon-index")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the deacon_indexes/ subdirectory will be created.",
)
@click.option(
    "--index-name",
    type=click.Choice(["panhuman-1", "panmouse-1"]),
    default="panhuman-1",
    show_default=True,
    help="Deacon index name to download.",
)
def get_deacon_index(path, index_name):
    """Download a pre-built Deacon minimizer index.

    Creates {path}/deacon_indexes/ and downloads the index using 'deacon index fetch'.
    """
    db_dir = Path(path) / "deacon_indexes"
    db_dir.mkdir(parents=True, exist_ok=True)

    output_file = db_dir / f"{index_name}.idx"

    cmd = ["deacon", "index", "fetch", index_name, "-o", str(output_file)]

    click.echo(f"Downloading Deacon index '{index_name}'...")
    _run(cmd)

    click.echo("\nDeacon index downloaded successfully.")
    click.echo(f"  Index: {output_file}")
    click.echo(f"Use --deacon-index {output_file} in your viralunity meta commands.")


@get_databases.command("all")
@click.option(
    "--path",
    default="databases",
    show_default=True,
    help="Parent directory where the databases will be created.",
)
@click.option(
    "--threads",
    default=4,
    show_default=True,
    type=int,
    help="Threads for diamond makedb.",
)
@click.option(
    "--refseq/--no-refseq",
    default=True,
    show_default=True,
    help="Limit to RefSeq genomes only.",
)
@click.pass_context
def get_all(ctx, path, threads, refseq):
    """Download kraken2, krona, taxdump and diamond databases at once."""
    click.echo(f"Downloading all databases to {path}...")
    ctx.invoke(
        get_kraken2,
        path=path,
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz",
    )
    ctx.invoke(get_krona, path=path)
    ctx.invoke(
        get_taxdump,
        path=path,
        url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    )
    ctx.invoke(
        get_diamond,
        path=path,
        taxon="Viruses",
        refseq=refseq,
        threads=threads,
        skip_makedb=False,
    )
    click.echo("\nAll databases downloaded successfully.")
