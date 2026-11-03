# ONT v3.3.0 Correctness Fixes

Date: 2026-11-03 (Asia/Shanghai).
Branch: codex/1028-ont-v330-performance.
Base: 6d28729ad9405e73a3930bbc356a74c5e97c1473.
Reason: fix the three confirmed bugs from the read-only branch review.

## Scope

Preserve ONT MAPQ/coverage thresholds, candidate ordering, primary alignment
selection, transcript IDs, workflow options and dependencies. Do not merge
unrelated legacy-backend changes from the 1103 branch.

1. Resolve missing transcript gene names from matching gene records, then exon
   records, then gene_id. Preserve explicit transcript names; do not inherit
   another transcript's alias. Reject a
   transcript without any usable gene identifier rather than silently writing
   an unassigned gene. Keep transcript/mRNA records as the candidate contract.
2. Search shifted segment offsets in indexed UMI candidate generation so the
   optimization retains Levenshtein neighbors with insertions/deletions. Keep
   directional abundance rules and exact bounded-distance verification.
3. Apply zero-total-cell removal even when feature filtering leaves zero rows.
   Preserve simultaneous initial cell/feature masks and normalization rules.

## Checklist

- [x] Add regressions and confirm failures against the unchanged branch.
- [x] Implement the three narrowly scoped fixes.
- [x] Run focused and full available Python/Rust tests, formatting and shell checks.
- [x] Review changes and record verification results and limitations here.

## Verification Record

Verified on 2026-11-03, Asia/Shanghai. The original failure cases were observed
before implementation, including an end-to-end synthetic gene-assignment TSV,
the 64/65-UMI indexing boundary and zero-feature matrix serialization.

- `cargo test --locked --offline --all-targets --quiet`: 108 passed.
- `python -m unittest discover -s tests -q`: 153 run, 150 passed, 3 skipped.
  Frozen BAM comparison used pysam 0.24.0's real samtools 1.23.1 dispatcher
  through a temporary command wrapper, not the Linux Conda samtools binary.
- `cargo build --locked --offline --bins`: passed; `flora --version` and
  `flora --help` work with unchanged version/options.
- `cargo fmt --all -- --check`, `git diff --check`, and `bash -n` for both
  workflows and both release packaging scripts passed.
- Added UMI distance and directional-abundance rejection tests in addition to
  exhaustive candidate completeness checks. Indexed candidate generation still
  verifies exact bounded Levenshtein distance before directional merging.
- Independent review found a fallback ordering ambiguity in the initial fix;
  a failing regression was added, explicit gene-over-exon precedence implemented,
  and the full suites rerun. Follow-up review found no remaining blockers.

The skipped tests require the external ONT assembly toolchain, original ONT
v3.3.0 oracle source, or Scanpy. Full assembly/remapping and large real Linux
datasets were not rerun here. Existing normal gene names retain their previous
labels; missing-name references may now use stable gene IDs in output tables.
No dependencies, default thresholds, version numbers or CLI options changed.

## History and Synchronization

Implementation is isolated in an independent workspace copy of the existing
branch under `Flora_private/.worktrees/1103-ont-v330-bugfixes`; no new branch
name is created. Record this operation with Git there for local rollback.
Synchronize only the five code/test files and this record to the same branch
in `/Users/liyy/Documents/GitHub/flora_private`. Leave that GitHub-local copy
uncommitted for the user's GitHub Desktop review; do not push or modify main.
The older dirty worktree is preserved untouched.
