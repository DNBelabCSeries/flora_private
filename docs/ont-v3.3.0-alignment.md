# ONT 3.3.0 comparison backend

Date: 2026-10-22 (Asia/Shanghai).
Branch: `codex/1022-ont-v330-feature-alignment`.
Base: `codex/1013-ont-prealignment-trimming`, commit `b931d16`.

## Reference identity

ONT withdrew wf-single-cell v3.3.0 because it contained critical issues.
The tag and release archive are unavailable. This branch uses the historical
source snapshot `d46f92121c79c28a02105cf2444029e91dc4b609` (2025-06-10),
where nextflow.config explicitly declares version `3.3.0` and CHANGELOG.md
starts with the 3.3.0 changes. It does not claim that this is the exact commit
of the deleted release tag. Source file SHA256 hashes and the original container
tag are recorded in `ont-v330-reference.json`.

Sources:

- [Historical source](https://github.com/epi2me-labs/wf-single-cell/tree/d46f92121c79c28a02105cf2444029e91dc4b609).
- [Current changelog with withdrawal notice](https://github.com/epi2me-labs/wf-single-cell/blob/v3.3.4/CHANGELOG.md).

The v3.3.4 comparison branch remains separate. Neither branch changes main.
No commits or pushes are made automatically; review and publish through
GitHub Desktop. The 3.3.0 development worktree is named
`1022-ont-v330-feature-alignment`; its files can be synchronized to the same
branch in the GitHub checkout without changing the separate 3.3.4 branch.

## Processing model

1. Retain Flora Glycine, barcode/UMI extraction, cell calls and cDNA FIFO.
2. Genome mapping: splice, forward-strand profile, MD tags, junction hints,
   10M batches and allocation cap; remove secondary records after alignment.
3. Assemble all mapped cDNA alignments per annotated chromosome, using
   StringTie 2.2.2 with `-L -c 2`, reference guidance and unique chromosome labels.
4. Extract transcript sequences using gffread; annotate with gffcompare.
5. Remap primary cDNA sequences to that chromosome's sample transcriptome with
   `map-ont --cap-kalloc 100m --cap-sw-mem 50m --end-bonus 10 -p 0.9 -N 3`.
   Keep transcriptome secondary alignments for candidate comparison.
6. Join candidate transcripts to reference annotations before unique/multi-hit
   classification. Unmatched novel transcripts remain in assembly exports but
   are not counted as reference features. Classes i/y/p/s lose transcript labels,
   but their gene labels remain eligible.
7. Unique surviving candidates require transcript MAPQ >0. Multiple candidates
   are sorted by alignment score, read coverage and transcript coverage. Distinct
   score/read coverage requires both coverages >=0.4; an otherwise tied candidate
   requires distinct transcript coverage >0.8. Gene assignment requires genomic
   MAPQ >=30; transcript assignment is not directly gated by genomic MAPQ.
8. Run Flora directional UMI correction using the new gene labels, then produce
   gene and transcript matrices. No per-cell/gene read cap is applied by default;
   `--cell-gene-max-reads` can impose one.

The default backend is `--feature-backend ont`. For the original Flora GTF
overlap backend use `--feature-backend legacy` and a separate output directory.
`--skip-isoform` skips transcript matrices, not assembly required for ONT gene
assignment. `--transcript-assign-mapq` applies to the legacy backend only.

## Comparison with 3.3.4

The inspected `assign_features.py` differs only in TSV compression. Genome and
transcriptome alignment options, coverage thresholds and candidate ranking are
the same. Matrix filtering, normalization to 10,000 and log10(1+x) are also the
same for non-spatial data. Thus this branch deliberately reuses the independently
implemented, tested selection engine instead of inventing version differences.

Differences in the reference sources include:

- 3.3.0 intermediate tables are plain TSV; 3.3.4 compresses several as zstd.
- 3.3.1 changed UMI-length filtering and sequencing saturation calculations.
- 3.3.1 updated StringTie to 2.2.3 to address assembly stalls; this historical
  comparison pins 2.2.2 and verifies its version before workflow execution.
- Later releases repaired several Visium HD and report edge cases.

Only the known StringTie version is pinned to the historical value. Other tools
use the existing Flora environment ranges, and actual versions are recorded.
The old ONT container is not executed, and no exact container/toolchain identity
or final numerical equivalence is claimed.

## Flora-specific boundaries

Flora's dual-end UMI representation can have lengths unlike a 10x UMI. Therefore
this branch does not reproduce the old ONT fixed-length UMI rejection behavior.
Flora's existing directional UMI engine, QC, sequencing saturation, clustering
and report definitions are retained. Visium, SNV and fusion workflows are not
added. Gene and transcript counting are independent, preserving transcript-only
reads rather than reproducing upstream in-place gene-dataframe filtering.

Count only primary mapped records; keep `(chromosome, read_id)` assignment keys
and reject duplicate primary read names. The original genomic BAM retains
supplementary records. Missing genes become positional UMI groups but are
excluded from the gene expression matrix, including prefixed mixed references.

The original TSV matrices remain. Added raw/processed MEX matrices contain
integer counts and normalized log values, respectively. Processed filtering uses
200 features/cell, 3 cells/feature, 20% mitochondrial counts and both MT-/mt-
identifiers for human/mouse data. Reports retain Flora cell calls, independently
of additional processed-MEX cell filtering. Empty processed results are exported
as empty matrices. These are deliberate integration choices, not a byte-exact
reproduction of ONT reports or preprocessing defaults.

## Output and failure handling

`transcriptome/` retains assembled FASTA/GTF, `tool_versions.json` (including
the historical commit and withdrawal flag), and `assignment_summary.json`.
Parameters TSV also records the reference commit. Assembly is sequential per
chromosome; SQLite read indexes avoid a whole-sample in-memory assignment map.
Temporary transcriptome BAMs are removed unless `--save-intermediate` or
`--full-output` is requested. Saved data use unique per-run directories.
Missing tools, wrong StringTie version and external failures stop explicitly.
No automatic fallback to the old annotation backend is performed.

## Test commands

The normal regression and real-tool smoke tests:

```bash
python -m unittest tests.test_ont_features tests.test_ont_matrices tests.test_ont_v330_reference -v
python -m unittest tests.test_ont_external -v
cargo test --locked --all-targets
```

Optional original-source oracle (research validation; source is not bundled):

```bash
mkdir -p /tmp/ont-v330-reference
curl -fL https://api.github.com/repos/epi2me-labs/wf-single-cell/tarball/d46f92121c79c28a02105cf2444029e91dc4b609 \
  -o /tmp/ont-v330-reference.tar.gz
tar -xzf /tmp/ont-v330-reference.tar.gz -C /tmp/ont-v330-reference --strip-components=1
python -m pip install polars pyarrow
ONT_V330_REFERENCE_DIR=/tmp/ont-v330-reference \
  python -m unittest tests.test_ont_v330_reference.OriginalSourceParityTests -v
```

The oracle checks all reference source hashes and executes the original ONT
parsing, joining, ranking and selection code on 500 generated reads. Only its
logging/parser helpers are supplied locally. Coverage, MAPQ, excluded classes,
novel candidates and read-group chunk boundaries vary in the generated data.
Exact-tie decisions also have focused unit tests in Flora. This verifies
feature selection, not end-to-end equivalence on the user's data.

## Implementation record

- [x] Recover and hash a verified version-3.3.0 historical source snapshot.
- [x] Compare assembly, mapping, assignment, matrix and UMI source changes.
- [x] Reuse verified independent engine; add reference identity and StringTie pin.
- [x] Original-source 500-read assignment oracle passes.
- [x] Run final regression, real-tool smoke, embedded build and code review.

Validation recorded on 2026-10-22: 74 Python checks passed, including the
500-read original-source oracle and real StringTie/gffread/gffcompare/minimap2
smoke test with repeated runs and the Rust UMI/matrix stages. All Rust targets,
shell syntax, formatting and whitespace checks passed. Release compilation and
the embedded single/mixed help passed. Local test versions: StringTie 2.2.2,
gffread 0.12.7, gffcompare 0.12.6 and minimap2 2.31-r1302; samtools 1.24 ran via
pysam's dispatcher. These are macOS synthetic checks; large Linux datasets,
full reports and full ONT Nextflow execution remain untested.
Code review also found and fixed unconditional ONT provenance in legacy-mode
parameter tables; both entrypoints now record the actual backend's identity.
