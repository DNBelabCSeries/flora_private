[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/Language-中文-0f766e)](README_zh-CN.md)

# Flora Private Development Repository

Private source, tests, and packaging tools for Flora. Keep this repository
private. The public `brilliantlee2/Flora` repository contains user documentation;
binary archives are distributed through GitHub Releases.

## Clone the private repository

```bash
git clone --branch codex/1028-ont-v330-performance --single-branch \
  https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

Access to the private repository is required. Use a GitHub token or authorized
SSH key; the SSH URL is `git@github.com:brilliantlee2/flora_private.git`.
Patched dependencies are included in `vendor/`; no submodule setup is needed.

## Build environment

Run these commands from the repository root on Linux x86_64:

```bash
conda env create -n flora-ont330 -f environment.yml
conda activate flora-ont330

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1
```

For an existing environment, use `conda env update -n flora-ont330 -f environment.yml`
instead of creating it, then reactivate it. Use Python 3.11.x.
`environment.yml` includes the build tools and runtime dependencies;
`environment.runtime.yml` is for packaged executables, and `requirements.txt`
lists Python packages only. This branch also requires StringTie 2.2.2,
`gffread`, and `gffcompare`, already included in both environment files. Keep
`flora-ont330` separate from the environment used for main or other ONT versions.

Before compiling, locate libclang for `rust-htslib`/`bindgen`. This is only a
build dependency, not a requirement for running the packaged executable.

```bash
LIBCLANG_SO="$(find -L "$CONDA_PREFIX/lib" -maxdepth 1 -type f \
  \( -name 'libclang.so' -o -name 'libclang.so.*' \) -print -quit)"
if [[ -n "$LIBCLANG_SO" ]]; then
  export LIBCLANG_PATH="$(dirname "$LIBCLANG_SO")"
  echo "LIBCLANG_PATH=$LIBCLANG_PATH"
else
  unset LIBCLANG_PATH
  echo "libclang not found. Update the Conda environment before building." >&2
fi
```

If not found, run `conda install -c conda-forge libclang`, then repeat the
lookup. Do not continue building until it succeeds. This example does not exit
your interactive SSH shell. Repeat the setup after opening a new shell, or
persist a successful value with
`conda env config vars set LIBCLANG_PATH="$LIBCLANG_PATH"` and reactivate `flora-ont330`.

## Development build and tests

Run the following checks once for each source revision you intend to package:

```bash
cargo test --locked
bash -n run_all.sh
bash -n run_all_mixed_species.sh
bash -n packaging/build_binary_release.sh
bash -n packaging/refresh_binary_release_metadata.sh
python -m unittest discover -s tests -v
```

Check that every command succeeds before continuing. Each `bash -n` checks one
script; passing multiple script paths to one invocation does not check them all.
These tests do not replace a full workflow run on representative data.

For a development executable without creating a release archive:

```bash
cargo build --release --locked
./target/release/flora --version
./target/release/flora --help
./target/release/flora mixed --help
./target/release/flora glycine --help
./target/release/flora analyze --help
```

Use `flora` for the single-species workflow, `flora mixed` for mixed species,
`flora glycine` for full-length identification only, and `flora analyze` for
barcode/UMI extraction and cell assignment only. Full workflows accept multiple
files after `--fastq`, or a non-recursive `--fastq-dir`. Glycine defaults to
10 concurrent jobs sharing 64 total threads.

## Build the Linux x86_64 release

After preparing the environment and passing the tests above, run:

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

The script performs its own locked release build; a separate `cargo build` is
not needed for packaging. It checks the Linux x86_64 binary and CLI, stages the
runtime environment and public documentation, and generates and verifies SHA256.
Python assets are embedded in `flora`, not shipped as loose source/bytecode files.

Outputs, with the version read from `Cargo.toml`:

```text
dist/Flora-<version>-linux-x86_64/
dist/Flora-<version>-linux-x86_64.tar.gz
dist/Flora-<version>-linux-x86_64.tar.gz.sha256
```

`BUILD_INFO.txt` inside the release records the required glibc version. Build
on a compatible Linux host; macOS cannot produce this release with this script.
An optional `FLORA_MAX_GLIBC=2.34` prefix rejects builds requiring newer symbols;
it does not lower the binary's glibc requirement. For older hosts, use a
compatible build environment or [Singularity](docs/singularity-build.md).

Before publishing, extract the archive into a fresh directory, activate the
runtime environment, and verify the executable independently:

```bash
./flora --version
./flora --help
./flora mixed --help
./flora glycine --help
./flora analyze --help
cat BUILD_INFO.txt
```

Also run a representative dataset before release. Help commands verify the CLI,
not the correctness of analysis results. The archive should contain only one
executable, `flora`, and no private source files.

## ONT comparison backend

This branch defaults to `--feature-backend ont`, an independent implementation
based on the historical [wf-single-cell v3.3.0 source snapshot](https://github.com/epi2me-labs/wf-single-cell/tree/d46f92121c79c28a02105cf2444029e91dc4b609).
ONT withdrew that version; this branch pins the source commit and StringTie
2.2.2 for comparison. It does not run ONT Nextflow or promise identical results.

The backend assembles a reference-guided sample transcriptome per chromosome,
remaps reads, then selects genes/transcripts by alignment score and coverage.
Gene assignment defaults to genomic MAPQ >=30; transcripts use separate
transcriptome evidence. Novel assemblies are exported, but only reference-linked
transcripts are counted. Flora retains its barcode/UMI model, cell calls, QC and
report definitions, without ONT's old 10x-specific fixed UMI-length filter.

Use the usual Flora command; the following options are optional:

| Option | Effect |
| --- | --- |
| `--feature-backend legacy` | Compare with the original GTF-overlap backend in a different output directory. |
| `--ont-temp-dir /local_scratch/flora_tmp` | Place disposable scratch files on a writable local disk; otherwise use the system temp directory/`TMPDIR`. |
| `--ont-jobs 2` | Run chromosomes concurrently; default: 1. |
| `--stringtie-opts '-c 2'` | Set StringTie assembly options. |
| `--cell-gene-max-reads 0` | Keep the ONT default of no per-cell/gene read truncation. |

Effective chromosome concurrency is capped at
`min(ont-jobs, max(1, cluster-threads / 8))`; jobs share the thread budget, with
additional tool/I/O helper threads possible. More jobs require more RAM and
scratch space. Avoid a small or memory-backed `/tmp` for large samples.
`--transcript-assign-mapq` and `--gene-tagging-backend` affect legacy only.
`--skip-isoform` skips transcript matrices, not the assembly needed for ONT genes.

In addition to the standard outputs:

- `transcriptome/<sample>.transcriptome.fa.gz` and `.transcriptome.gtf.gz`:
  assembled transcript sequences and annotation.
- `transcriptome/tool_versions.json` and `assignment_summary.json`:
  actual tool versions and assignment-status counts.
- Gene/transcript raw and processed MEX directories under `matrix/`, plus
  `<sample>.matrix_processing.json` with processing details.

Raw MEX contains integer UMI counts. Processed MEX applies minimum 200 features
per cell, minimum 3 cells per feature, and maximum 20% mitochondrial counts,
then normalizes to 10,000 counts/cell and applies `log10(1+x)`. Both `MT-` and
`mt-` identifiers are recognized. The report uses Flora's original cell calls,
not the additionally filtered processed-MEX population. Temporary transcriptome
BAMs are retained only with `--save-intermediate` or `--full-output`.

Implementation and performance details:
[ONT backend](docs/ont-v3.3.0-alignment.md),
[performance record](docs/ont-v330-performance.md).

## Workflow output policy

Default light output retains the final tagged BAM/index, expression matrices,
cell/barcode mapping, RNA clustering, QC, `qc/metrics_summary.xlsx`, logs, and
the HTML report. After successful completion, reproducible read-level tables,
the pre-tagging aligned BAM, and assignment TSVs are removed. An empty
`alignment/` directory is also removed.

| Option | Effect |
| --- | --- |
| `--save-intermediate` | Retain read-level assignments and other debugging intermediates. |
| `--full-output` | Retain upstream FASTQ outputs and workflow intermediates. |
| `--remove-final-bam` | Remove all BAM/index files in `alignment/` and `matrix/` after successful analysis, including the final tagged BAM. |

These switches control retained files, not analytical thresholds. BAM removal
takes precedence when combined with either retention switch. Full workflows do
not enable Glycine's extra classification FASTQs; run `flora glycine` with
`--keep-all-outputs` to obtain those files.

## Publish a release

1. Update `Cargo.toml`, then run `cargo check` without `--locked` to update the
   Flora package version in `Cargo.lock`.
2. Update versioned public documentation, release notes, and
   `packaging/singularity/Flora.def` together.
3. Run tests, build the release, and validate it as described above.
4. Sync the approved files under `docs/repository-templates/public/` to the
   public repository. Upload `.tar.gz` and `.sha256` as Release assets; an SIF
   and its checksum are optional assets.

Use the same version for the tag and title, for example `v0.1.2` and
`Flora v0.1.2`. Do not commit archives, SIF images, Docker tar files, Miniforge
installers, or `dist/` to Git. Never copy private `src/`, `scripts/`, `vendor/`,
Cargo files, tests, or packaging code to the public repository. A compiled
executable is not encryption and may still be reverse engineered.

## Source layout

| Path | Purpose |
| --- | --- |
| `src/`, `vendor/`, `Cargo.toml`, `Cargo.lock` | Rust source and pinned dependencies. |
| `scripts/` | Python analysis, QC, and reporting assets. |
| `run_all.sh`, `run_all_mixed_species.sh` | Single- and mixed-species workflow scripts. |
| `tests/` | Regression and release checks. |
| `packaging/` | Binary and container packaging tools. |
| `docs/repository-templates/` | Public/private repository document templates. |

Detailed container instructions: [Build a Singularity image](docs/singularity-build.md).
