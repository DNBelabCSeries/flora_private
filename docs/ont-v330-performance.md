# ONT 3.3.0 Performance Plan and Operation Record

Date: 2026-10-28 (Asia/Shanghai).
Base: codex/1022-ont-v330-feature-alignment, d4d0fe6.
Branch: codex/1028-ont-v330-performance.

## Contract

Keep feature selection, MAPQ/coverage thresholds, candidate tie order, duplicate
primary-ID validation, primary/supplementary filtering, TSV order/columns,
GN/TX tags, UMI clustering, matrices and reporting unchanged. Do not transplant
ONT's fixed UMI length filter. BAM compression bytes are not an equivalence
criterion; decoded records and tags are.

## First-Round Implementation

- [x] Add parity and performance-structure regression tests before implementation.
- [x] Batch SQLite joins and duplicate-ID validation with bounded record batches.
      Use disposable, non-resumable scratch databases without durability fsyncs;
      final outputs are not used as databases and retain normal file semantics.
- [x] Replace per-base Python quality conversion with pysam's native conversion;
      buffer FASTQ writes and preserve missing-quality and reverse-strand behavior.
- [x] Add --ont-temp-dir (default: system temp/TMPDIR) for scratch storage.
- [x] Add --ont-jobs (default: 1), cap effective jobs by total --cluster-threads,
      and use an isolated database/directory per chromosome. Merge in original
      annotation order, never completion order, with bounded scheduling.
- [x] Log feature joins, final output, indexing, chromosome completion and rates.
- [x] Wire both workflows, embedded CLI options, fixtures and bilingual READMEs.
- [x] Run old/new decoded BAM and byte-exact TSV comparisons, real-tool smoke
      tests where available, Rust build and focused workflow/release tests.

Removing intermediate CB/UR/GN BAMs is deferred. The approved second round below
adds upstream coordinate reuse without changing the extracted barcode fields.

## Verification

No commit or push is performed automatically; the user reviews and publishes
through GitHub Desktop.

- 84 focused Python tests passed, including the 500-read original ONT source
  oracle and hashes, real assembly/remapping, serial-versus-two-job comparison,
  saved-intermediate reruns, and Rust UMI/matrix integration.
- All Rust targets passed; release executable builds and exposes both options.
- Real smoke tools: StringTie 2.2.2, minimap2 2.31-r1302, gffread 0.12.7,
  gffcompare 0.12.6. Samtools commands in this local macOS test used pysam's
  embedded samtools 1.23.1 dispatcher, not the Linux Conda-pinned 1.24 binary.
- Independent review found no remaining blockers; an additional 20,031-read
  batch-boundary comparison matched TSVs, decoded BAM tags and summary counts.
- A one-million-read, 1,000-base synthetic benchmark against d4d0fe6 matched
  both TSVs byte-for-byte and all decoded BAM records. FASTQ conversion alone:
  30.01s before, 2.27s after (13.24x). Final BAM/TSV output: 9.07s before,
  9.15s after on warm local storage, i.e. no demonstrated output speedup there.
  This does NOT predict end-to-end speedup or shared-storage performance.
- The scratch database uses MEMORY journaling and synchronous=OFF. MEMORY
  journaling retains SAVEPOINT rollback for correct duplicate-ID diagnostics;
  journal_mode=OFF is explicitly avoided. These databases are disposable and
  must not be reused after a process/host failure.

Large real Linux datasets and peak memory under concurrent chromosome work
still require server-side validation. Defaults keep one chromosome job.

## Second Round: Same Branch

Approved continuation on 2026-10-28; no additional branch or automatic commit.

- [x] Cache two independent terminal bounds per FASTQ record in an anonymous
  temporary binary file (16 bytes/read, not a sample-sized RAM array). Apply the
  final barcode-validity decisions only when writing cDNA; retain the old
  coordinate search for the first few fallback diagnostics and parity tests.
- [x] Compute Read QC in bounded parallel chunks and merge in input order.
- [x] Avoid intermediate trimmed-record String copies and log input versus
  output/pipe-backpressure time separately.
- [x] Retain isolated chromosome feature databases instead of copying every
  row to a global feature database. Keep at most one read-only shard open;
  global duplicate-ID detection stays separate and unchanged.
- [x] Delete chromosome MAPQ scratch separately after assignment. Preserve
  output order, saved-intermediate handling and error behavior.
- [x] Compare cached/uncached cDNA and parallel/serial QC exactly; validate
  old/new database layouts, BAM tags, TSV bytes and matrices; rebuild executable.

### Second-Round Verification

Operation recorded: 2026-10-28 20:39 CST. Reason: continue reducing repeated
terminal searches and disposable database I/O on the same branch.

- 103 Rust tests across all targets and 87 focused Python tests passed. The
  external tool/original-source caveats listed above still apply.
- Cached bounds match the original interval function for all four corrected-end
  validity combinations, including negative/out-of-range coordinates. Cached
  and uncached FASTQ output is byte-exact with input gaps, duplicate physical
  IDs, different batch boundaries and mixed gzip/plain input files.
- Sequential and parallel Read QC JSON match exactly with varying read lengths
  and quality strings, including across chunk/batch boundaries.
- A 30,000-read, 1,082-base, four-thread release-executable comparison against
  the first-round binary matched 20 upstream text/JSON/cDNA files byte-for-byte.
  Its sub-second local timings are too short to extrapolate to real datasets.
- A 20,006-record writer comparison crosses both batch and chromosome boundaries:
  global-index and per-chromosome lookup give byte-exact TSVs, identical decoded
  BAM records/tags and identical summary counts. Global duplicate checks remain.
- Separate MAPQ scratch, read-only feature connections, worker scratch cleanup,
  saved intermediates, release build and shell syntax were checked.
- Review identified two allocation regressions, then failing regression tests
  confirmed both: reserve the full input batch before appending 512-row chunks
  (preserves the previous long-lived row-vector capacity), and split the worker's
  32 MiB SQLite cache target into 16 MiB per schema instead of doubling it.

The coordinate cache needs 16 bytes per input read on temporary storage (about
1.49 GiB per 100 million reads), not an additional all-read RAM array. Prefer
a local disk-backed SSD via `--ont-temp-dir`; a tmpfs directory charges this space
to memory. It is closed/deleted immediately after cDNA output (also on errors).

Intermediate BAM fusion and elimination of the second FASTQ decompression pass
remain out of scope; the second pass now avoids repeated terminal searches.
