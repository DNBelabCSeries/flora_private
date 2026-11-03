"""Flora feature assignment using wf-single-cell v3.3.0 historical rules.

Independent implementation. Keep read-sized indexes on disk and process one
chromosome at a time; transcriptome secondary alignments are needed for ranking.
"""

from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
import gzip
import itertools
import json
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time


REFERENCE_VERSION = "wf-single-cell v3.3.0"
REFERENCE_COMMIT = "d46f92121c79c28a02105cf2444029e91dc4b609"
REFERENCE_CONTAINER = "ontresearch/wf-single-cell:shab5b0dea0efc4685f74c8b4f91c979c587e23a020"
REFERENCE_WITHDRAWN = True
EXCLUDED_CLASSES = {"i", "y", "p", "s"}
BATCH_SIZE = 10000


@dataclass(frozen=True)
class Hit:
    gene: str
    transcript: str
    score: int
    query_coverage: float
    transcript_coverage: float
    mapq: int


def select_features(hits, genome_mapq, min_mapq=30, min_read_cov=0.4, min_tr_cov=0.4):
    if not hits:
        return "-", "-"
    ranked = sorted(hits, key=lambda h: (h.score, h.query_coverage, h.transcript_coverage), reverse=True)
    best = ranked[0]
    gene = best.gene if genome_mapq >= min_mapq else "-"
    if len(ranked) == 1:
        return (gene, best.transcript) if best.mapq > 0 else ("-", "-")
    second = ranked[1]
    score_or_query_diff = best.score != second.score or best.query_coverage != second.query_coverage
    sufficient_coverage = best.query_coverage >= min_read_cov and best.transcript_coverage >= min_tr_cov
    coverage_tiebreak = best.transcript_coverage != second.transcript_coverage and best.transcript_coverage > 0.8
    keep = sufficient_coverage if score_or_query_diff else coverage_tiebreak
    return gene, best.transcript if keep else "-"


def log(message):
    print(f"[ont_features] {message}", file=sys.stderr, flush=True)


def open_scratch_database(path, genome_path=None):
    db = sqlite3.connect(path)
    # Disposable scratch only: a failed run is rebuilt, never resumed from this DB.
    db.execute("PRAGMA journal_mode=MEMORY")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA cache_size=-32768")
    db.execute("PRAGMA temp_store=MEMORY")
    db.execute("CREATE TABLE features(chrom TEXT,read_id TEXT,gene TEXT,transcript TEXT,reason TEXT,PRIMARY KEY(chrom,read_id)) WITHOUT ROWID")
    if genome_path is None:
        db.execute("CREATE TABLE genome(read_id TEXT PRIMARY KEY,mapq INTEGER) WITHOUT ROWID")
    else:
        # MAPQ lookup data can be discarded without retaining/freeing its pages in features.
        db.execute("ATTACH DATABASE ? AS maps", (str(genome_path),))
        db.execute("PRAGMA maps.journal_mode=MEMORY")
        db.execute("PRAGMA maps.synchronous=OFF")
        db.execute("PRAGMA main.cache_size=-16384")
        db.execute("PRAGMA maps.cache_size=-16384")
        db.execute("CREATE TABLE maps.genome(read_id TEXT PRIMARY KEY,mapq INTEGER) WITHOUT ROWID")
    return db


class ChromosomeFeatureLookup:
    """Keep at most one read-only chromosome index open during coordinate-order output."""

    def __init__(self, paths):
        self.paths = paths
        self.chrom = None
        self.db = None

    def lookup(self, chrom, names):
        path = self.paths.get(chrom)
        if path is None:
            return {}
        if chrom != self.chrom:
            self.close()
            self.db = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)
            self.db.execute("PRAGMA cache_size=-32768")
            self.db.execute("PRAGMA temp_store=MEMORY")
            self.chrom = chrom
        return lookup_batch(self.db, "features", names, chrom)

    def close(self):
        if self.db is not None:
            self.db.close()
        self.chrom = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def lookup_batch(db, table, names, chrom=None):
    if table not in {"features", "genome"}:
        raise ValueError(f"Invalid lookup table: {table}")
    db.execute("CREATE TEMP TABLE IF NOT EXISTS lookup_reads(read_id TEXT PRIMARY KEY) WITHOUT ROWID")
    db.execute("DELETE FROM lookup_reads")
    db.executemany("INSERT INTO lookup_reads VALUES (?)", ((n,) for n in sorted(set(names))))
    # Read requests in key order for locality; one join per batch, not one SELECT/read.
    if table == "genome":
        rows = db.execute("SELECT r.read_id,g.mapq FROM lookup_reads r CROSS JOIN genome g ON g.read_id=r.read_id")
    else:
        rows = db.execute("SELECT r.read_id,f.gene,f.transcript,f.reason FROM lookup_reads r CROSS JOIN features f ON f.read_id=r.read_id AND f.chrom=?", (chrom,))
    return {row[0]: row[1:] for row in rows}


def validate_read_ids(db, names):
    db.execute("SAVEPOINT validate_ids")
    try:
        db.executemany("INSERT INTO counted_reads VALUES (?)", ((n,) for n in sorted(names)))
    except sqlite3.IntegrityError:
        db.execute("ROLLBACK TO validate_ids")
        # Error path only: replay to report the same offending input name.
        for name in names:
            try:
                db.execute("INSERT INTO counted_reads VALUES (?)", (name,))
            except sqlite3.IntegrityError as error:
                raise ValueError(f"Duplicate primary read name: {name}; input FASTQ read IDs must be unique") from error
        raise
    finally:
        db.execute("RELEASE validate_ids")


def ordered_parallel_map(function, items, jobs):
    if jobs == 1:
        yield from map(function, items)
        return
    items = iter(items)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        pending = deque(pool.submit(function, item) for item in itertools.islice(items, jobs))
        try:
            while pending:
                result = pending.popleft().result()
                item = next(items, None)
                if item is not None:
                    pending.append(pool.submit(function, item))
                yield result
        finally:
            for future in pending:
                future.cancel()


def fastq_bytes(record):
    import pysam

    name = record.query_name
    sequence = record.get_forward_sequence()
    if sequence is None:
        raise ValueError(f"Missing sequence for {name}")
    qualities = record.get_forward_qualities() if record.query_qualities is not None else None
    quality = pysam.array_to_qualitystring(qualities) if qualities is not None else "B" * len(sequence)
    return f"@{name}\n{sequence}\n+\n{quality}\n".encode()


def execution_layout(threads, jobs):
    if threads < 1 or jobs < 1:
        raise ValueError("threads and jobs must be positive")
    effective = min(jobs, max(1, threads // 8))
    return effective, threads // effective


def process_chromosome(item, args, root, genome, threads):
    number, chrom, gtf = item
    started = time.monotonic()
    directory = root / f"chrom_{number:05d}"
    directory.mkdir()
    genome_db = directory / "genome.sqlite"
    db = open_scratch_database(directory / "features.sqlite", genome_path=genome_db)
    try:
        assemble(args.bam, chrom, gtf, genome, directory, threads, args.stringtie_opts, f"flora_chr{number}")
        fasta = directory / "transcriptome.fa"
        if fasta.exists() and fasta.stat().st_size:
            remap_chromosome(db, args.bam, chrom, fasta, directory / "transcript.bam", threads)
            candidates = load_candidates(directory / "gffcompare.stringtie.gff.tmap", gtf)
            assign_chromosome(db, chrom, directory / "transcript.bam", candidates, args.mapq)
        db.commit()
    finally:
        db.close()
    genome_db.unlink()
    log(f"Stage chromosome ({chrom}) elapsed: {time.monotonic() - started:.2f}s (includes its internal stages)")
    return directory, chrom, gtf


def timed_run(command):
    started = time.monotonic()
    log(shlex.join([str(x) for x in command]))
    subprocess.run(command, check=True)
    log(f"Stage {command[0]} elapsed: {time.monotonic() - started:.2f}s")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attributes(text):
    return dict(re.findall(r'(\w+)\s+"([^"]*)"', text))


def split_annotation(gtf, root):
    paths = {}
    last_chrom = None
    handle = None
    try:
        with open_text(gtf) as source:
            for line in source:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError("Reference annotation must be a nine-column GTF")
                chrom = fields[0]
                if chrom not in paths:
                    paths[chrom] = root / f"chrom_{len(paths):05d}.gtf"
                if chrom != last_chrom:
                    if handle:
                        handle.close()
                    handle = paths[chrom].open("a")
                    last_chrom = chrom
                handle.write(line)
    finally:
        if handle:
            handle.close()
    return paths


def reference_transcripts(gtf):
    transcripts, gene_names, exon_names = {}, {}, {}
    with open_text(gtf) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            attr = attributes(fields[8])
            feature = fields[2].lower()
            gene_id = attr.get("gene_id")
            if gene_id in {"", "-", ".", "NA"}:
                gene_id = None
            name = next((attr[key] for key in ("gene_name", "genename")
                         if attr.get(key) not in {None, "", "-", ".", "NA"}), None)
            if gene_id and name:
                if feature == "gene":
                    gene_names.setdefault(gene_id, name)
                elif feature == "exon":
                    exon_names.setdefault(gene_id, name)
            tx = attr.get("transcript_id")
            if tx and feature in {"transcript", "mrna"}:
                transcripts[tx] = (gene_id, name)
    result = {}
    # Resolve after the full scan: gene/exon names may follow transcript records.
    for tx, (gene_id, name) in transcripts.items():
        gene = name or gene_names.get(gene_id) or exon_names.get(gene_id) or gene_id
        if not gene:
            raise ValueError(f"Reference transcript {tx} has no gene_name or gene_id in {gtf}")
        result[tx] = gene
    return result


def load_candidates(tmap, gtf):
    reference = reference_transcripts(gtf)
    result = {}
    with open(tmap) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx = row["ref_id"]
            # ONT counts only assembled candidates linked to reference transcripts.
            if tx not in reference:
                continue
            assigned_tx = "-" if row["class_code"] in EXCLUDED_CLASSES else tx
            result[row["qry_id"]] = (reference[tx], assigned_tx)
    return result


def make_hit(record, candidates, lengths):
    pair = candidates.get(record.reference_name)
    if pair is None:
        return None
    read_length = record.infer_read_length() or 0
    if not read_length:
        return None
    return Hit(
        pair[0], pair[1], record.get_tag("AS") if record.has_tag("AS") else 0,
        record.query_alignment_length / read_length,
        (record.reference_end - record.reference_start) / lengths[record.reference_name],
        record.mapping_quality,
    )


def assign_chromosome(db, chrom, transcript_bam, candidates, min_mapq):
    import pysam

    started = time.monotonic()
    processed = 0
    last_progress = started
    with pysam.AlignmentFile(transcript_bam, "rb", check_sq=False) as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        batch = []
        def flush():
            mapqs = lookup_batch(db, "genome", [name for name, _ in batch])
            pending = []
            for name, hits in batch:
                if name not in mapqs:
                    continue
                gene, tx = select_features(hits, mapqs[name][0], min_mapq)
                reason = "candidate" if hits else "no_reference_candidate"
                if hits and gene == "-" and tx == "-":
                    reason = "ambiguous_or_zero_transcript_mapq"
                pending.append((chrom, name, gene, tx, reason))
            db.executemany("INSERT OR REPLACE INTO features VALUES (?,?,?,?,?)", pending)
            db.commit()
            batch.clear()
        for name, records in itertools.groupby(bam, key=lambda a: a.query_name):
            hits = [make_hit(a, candidates, lengths) for a in records
                    if not a.is_unmapped and not a.is_supplementary]
            hits = [h for h in hits if h is not None]
            batch.append((name, hits))
            processed += 1
            if len(batch) >= BATCH_SIZE:
                flush()
                now = time.monotonic()
                if processed % 1000000 == 0 or now - last_progress >= 30:
                    log(f"feature_assignment ({chrom}): {processed} reads, {processed / max(now - started, 0.001):.0f} reads/s")
                    last_progress = now
        if batch:
            flush()
    log(f"Stage feature_assignment ({chrom}) elapsed: {time.monotonic() - started:.2f}s; reads={processed}")


def assemble(bam, chrom, gtf, genome, directory, threads, options, label):
    started = time.monotonic()
    view = subprocess.Popen(["samtools", "view", "-h", str(bam), chrom], stdout=subprocess.PIPE)
    command = ["stringtie", "-L", *shlex.split(options), "-p", str(threads),
               "-G", str(gtf), "-l", label, "-o", str(directory / "stringtie.gff"), "-"]
    log(shlex.join(command))
    try:
        result = subprocess.run(command, stdin=view.stdout)
    finally:
        view.stdout.close()
        view_code = view.wait()
    if view_code:
        raise subprocess.CalledProcessError(view_code, view.args)
    result.check_returncode()
    log(f"Stage stringtie ({chrom}) elapsed: {time.monotonic() - started:.2f}s")
    # An empty assembly is a valid outcome for sparse chromosomes.
    with (directory / "stringtie.gff").open() as handle:
        has_exons = any(not line.startswith("#") and "\texon\t" in line for line in handle)
    if not has_exons:
        log(f"No assembled transcripts on {chrom}")
        return
    timed_run(["gffread", "-g", str(genome), "-w", str(directory / "transcriptome.fa"),
               str(directory / "stringtie.gff")])
    timed_run(["gffcompare", "-o", str(directory / "gffcompare"), "-r", str(gtf),
               str(directory / "stringtie.gff")])


def remap_chromosome(db, bam_path, chrom, fasta, output, threads):
    import pysam

    started = time.monotonic()
    sort_threads = min(3, max(1, threads // 4))
    view_threads = int(threads >= 8)
    mm = ["minimap2", "-ax", "map-ont", "--cap-kalloc", "100m", "--cap-sw-mem", "50m",
          "--end-bonus", "10", "-p", "0.9", "-N", "3", "-t", str(max(threads - sort_threads - view_threads, 1)), str(fasta), "-"]
    log(shlex.join(mm))
    processes = []
    try:
        aligner = subprocess.Popen(mm, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        processes.append(aligner)
        viewer = subprocess.Popen(["samtools", "view", "-h", "-@", str(view_threads), "-u", "-F", "2052", "-"],
                                  stdin=aligner.stdout, stdout=subprocess.PIPE)
        processes.append(viewer)
        aligner.stdout.close()
        sorter = subprocess.Popen(["samtools", "sort", "-n", "-@", str(sort_threads), "--no-PG", "-o", str(output), "-"],
                                  stdin=viewer.stdout)
        processes.append(sorter)
        viewer.stdout.close()
        pending = []
        buffer = bytearray()
        n_reads = 0
        with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
            for record in bam.fetch(chrom):
                if record.is_unmapped or record.is_secondary or record.is_supplementary:
                    continue
                buffer.extend(fastq_bytes(record))
                if len(buffer) >= 1024 * 1024:
                    aligner.stdin.write(buffer)
                    buffer.clear()
                pending.append((record.query_name, record.mapping_quality))
                n_reads += 1
                if len(pending) >= BATCH_SIZE:
                    db.executemany("INSERT OR REPLACE INTO genome VALUES (?,?)", pending)
                    pending.clear()
        if buffer:
            aligner.stdin.write(buffer)
        db.executemany("INSERT OR REPLACE INTO genome VALUES (?,?)", pending)
        db.commit()
        aligner.stdin.close()
        log(f"Stage transcript_input ({chrom}) elapsed: {time.monotonic() - started:.2f}s; reads={n_reads} (includes pipe backpressure)")
        codes = [p.wait() for p in processes]
        for p, code in zip(processes, codes):
            if code:
                raise subprocess.CalledProcessError(code, p.args)
    except BaseException:
        for p in processes:
            if p.poll() is None:
                p.terminate()
        for p in processes:
            p.wait()
        raise
    log(f"Stage transcript_remapping ({chrom}) elapsed: {time.monotonic() - started:.2f}s")


def write_tagged_assignments(db, tagged_bam, output_bam, gene_tsv, transcript_tsv, min_mapq, threads,
                             feature_lookup=None):
    import pysam
    from collections import Counter

    stats = Counter()
    started = time.monotonic()
    processed = 0
    lookup_seconds = 0.0
    last_progress = started
    log("Starting tagged_assignment_output: bounded batches, primary records only")
    db.execute("CREATE TABLE IF NOT EXISTS counted_reads(read_id TEXT PRIMARY KEY) WITHOUT ROWID")
    with pysam.AlignmentFile(tagged_bam, "rb", threads=max(1, min(threads // 4, 2))) as source, \
            pysam.AlignmentFile(output_bam, "wb", template=source, threads=max(1, min(threads, 8))) as dest, \
            open(gene_tsv, "w") as genes, open(transcript_tsv, "w") as transcripts:
        primary = (r for r in source if not (r.is_secondary or r.is_supplementary or r.is_unmapped))
        while batch := list(itertools.islice(primary, BATCH_SIZE)):
            phase = time.monotonic()
            validate_read_ids(db, [r.query_name for r in batch])
            by_chrom = {}
            for record in batch:
                by_chrom.setdefault(record.reference_name, []).append(record.query_name)
            features = {chrom: (feature_lookup.lookup(chrom, names) if feature_lookup is not None
                                else lookup_batch(db, "features", names, chrom))
                        for chrom, names in by_chrom.items()}
            lookup_seconds += time.monotonic() - phase
            gene_lines, tx_lines = [], []
            for record in batch:
                rid = record.query_name
                gene, tx, reason = features[record.reference_name].get(rid, ("-", "-", "no_reference_candidate"))
                status = "Assigned" if gene != "-" else (
                    "Unassigned_mapq" if record.mapping_quality < min_mapq else "Unassigned_no_features")
                if gene == "-" and reason == "ambiguous_or_zero_transcript_mapq" and record.mapping_quality >= min_mapq:
                    status = "Unassigned_ambiguous"
                tx_status = "Assigned" if tx != "-" else (
                    "Unassigned_ambiguous" if reason != "no_reference_candidate" else "Unassigned_no_features")
                gene_out, tx_out = (gene if gene != "-" else "NA"), (tx if tx != "-" else "NA")
                gene_lines.append(f"{rid}\t{status}\t{record.mapping_quality}\t{gene_out}\n")
                tx_lines.append(f"{rid}\t{tx_status}\t{record.mapping_quality}\t{gene_out}\t{tx_out}\n")
                record.set_tag("GN", gene_out, value_type="Z")
                record.set_tag("TX", tx_out, value_type="Z")
                dest.write(record)
                stats[f"gene:{status}"] += 1
                stats[f"transcript:{tx_status}"] += 1
                stats[f"reason:{reason}"] += 1
            genes.writelines(gene_lines)
            transcripts.writelines(tx_lines)
            db.commit()
            processed += len(batch)
            now = time.monotonic()
            if processed % 1000000 == 0 or now - last_progress >= 30:
                log(f"tagged_assignment_output: {processed} reads, {processed / max(now - started, 0.001):.0f} reads/s")
                last_progress = now
    log(f"Stage batched_feature_lookup elapsed: {lookup_seconds:.2f}s (included in tagged_assignment_output)")
    log(f"Stage tagged_assignment_output elapsed: {time.monotonic() - started:.2f}s; reads={processed}")
    started = time.monotonic()
    pysam.index("-@", str(max(1, min(threads, 8))), str(output_bam))
    log(f"Stage feature_bam_index elapsed: {time.monotonic() - started:.2f}s")
    return dict(stats)


def run(args):
    import pysam

    for name in ("minimap2", "samtools", "stringtie", "gffread", "gffcompare"):
        if shutil.which(name) is None:
            raise RuntimeError(f"Missing required tool: {name}; update the Flora conda environment")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    versions = {"behavior_reference": REFERENCE_VERSION, "reference_commit": REFERENCE_COMMIT,
                "reference_container": REFERENCE_CONTAINER, "reference_withdrawn": REFERENCE_WITHDRAWN}
    for name, flag in [("minimap2", "--version"), ("samtools", "--version"),
                       ("stringtie", "--version"), ("gffread", "--version"), ("gffcompare", "--version")]:
        p = subprocess.run([name, flag], capture_output=True, text=True, errors="replace")
        versions[name] = ((p.stdout + p.stderr).strip().splitlines() or ["unknown"])[0]
    (out / "tool_versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    if versions["stringtie"] != "2.2.2":
        raise RuntimeError("The ONT 3.3.0 comparison backend requires StringTie 2.2.2; update the pinned Flora conda environment")
    log("ONT 3.3.0 was withdrawn upstream; using the recorded historical source commit for comparison")
    log(f"Behavior reference: {REFERENCE_VERSION}; {json.dumps(versions)}")
    jobs, worker_threads = execution_layout(args.threads, getattr(args, "jobs", 1))
    scratch_parent = getattr(args, "temp_dir", None)
    if scratch_parent is not None:
        scratch_parent = Path(scratch_parent).resolve()
        scratch_parent.mkdir(parents=True, exist_ok=True)
    log(f"Scheduling: requested_jobs={getattr(args, 'jobs', 1)} effective_jobs={jobs} total_threads={args.threads} threads_per_job={worker_threads}; auxiliary tool threads may be additional")
    with tempfile.TemporaryDirectory(prefix="ont-work-", dir=scratch_parent) as temp:
        root = Path(temp)
        log(f"Scratch directory: {root}; disposable databases are not resumable")
        free_gib = shutil.disk_usage(root).free / (1024 ** 3)
        log(f"Scratch free space: {free_gib:.1f} GiB; check capacity for large datasets")
        genome = args.genome.resolve()
        if not Path(str(genome) + ".fai").exists():
            # Each invocation gets its own link/index, including reruns with a new reference.
            linked_genome = root / "reference.fa"
            linked_genome.symlink_to(genome)
            timed_run(["samtools", "faidx", str(linked_genome)])
            genome = linked_genome
        saved_run = out / "intermediates" / root.name
        if args.save_intermediate:
            saved_run.mkdir(parents=True)
        db = open_scratch_database(root / "assignments.sqlite")
        try:
            feature_root = root / "feature_index"
            feature_root.mkdir()
            feature_paths = {}
            annotations = split_annotation(args.gtf, root)
            if not reference_transcripts(args.gtf):
                raise ValueError("ONT feature assignment requires reference transcript/mRNA records with transcript_id and gene_name or gene_id in the GTF")
            with pysam.AlignmentFile(args.bam, "rb") as bam:
                mapped = {s.contig for s in bam.get_index_statistics() if s.mapped}
            transcriptome_fa = out / f"{args.sample}.transcriptome.fa.gz"
            transcriptome_gtf = out / f"{args.sample}.transcriptome.gtf.gz"
            work = [(n, chrom, gtf) for n, (chrom, gtf) in enumerate(annotations.items()) if chrom in mapped]
            def worker(item):
                return process_chromosome(item, args, root, genome, worker_threads)
            with gzip.open(transcriptome_fa, "wt") as fa_out, gzip.open(transcriptome_gtf, "wt") as gtf_out:
                for directory, chrom, gtf in ordered_parallel_map(worker, work, jobs):
                    fasta = directory / "transcriptome.fa"
                    annotated = directory / "gffcompare.annotated.gtf"
                    if fasta.exists():
                        with fasta.open() as handle:
                            shutil.copyfileobj(handle, fa_out)
                    if annotated.exists():
                        with annotated.open() as handle:
                            shutil.copyfileobj(handle, gtf_out)
                    phase = time.monotonic()
                    feature_paths[chrom] = feature_root / f"{directory.name}.sqlite"
                    (directory / "features.sqlite").rename(feature_paths[chrom])
                    log(f"Stage register_feature_index ({chrom}) elapsed: {time.monotonic() - phase:.2f}s; no database copy")
                    if args.save_intermediate:
                        shutil.copyfile(gtf, directory / "reference.gtf")
                        shutil.move(str(directory), str(saved_run / directory.name))
                    else:
                        shutil.rmtree(directory)
            with ChromosomeFeatureLookup(feature_paths) as feature_lookup:
                stats = write_tagged_assignments(db, args.tagged_bam, args.output_bam, args.gene_assigns,
                                                 args.transcript_assigns, args.mapq, args.threads,
                                                 feature_lookup=feature_lookup)
            (out / "assignment_summary.json").write_text(json.dumps(stats, indent=2) + "\n")
        finally:
            db.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("bam", "tagged-bam", "genome", "gtf", "output-dir", "output-bam", "gene-assigns", "transcript-assigns"):
        p.add_argument(f"--{name}", required=True, type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--mapq", type=int, default=30)
    p.add_argument("--stringtie-opts", default="-c 2")
    p.add_argument("--save-intermediate", action="store_true")
    p.add_argument("--temp-dir", type=Path, help="Scratch parent; defaults to the system temp directory/TMPDIR")
    p.add_argument("--jobs", type=int, default=1, help="Concurrent chromosome jobs sharing --threads; default: 1")
    args = p.parse_args()
    if args.threads < 1 or args.jobs < 1 or not 0 <= args.mapq <= 255:
        p.error("threads/jobs must be positive and MAPQ must be between 0 and 255")
    run(args)


if __name__ == "__main__":
    main()
