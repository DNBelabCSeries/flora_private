[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/Language-中文-0f766e)](README_zh-CN.md)

# Flora Private Development Repository

Private source, tests, and packaging tools for Flora. Keep this repository
private. The public `brilliantlee2/Flora` repository contains user documentation;
binary archives are distributed through GitHub Releases.

## Clone the private repository

```bash
git clone https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

Access to the private repository is required. Use a GitHub token or authorized
SSH key; the SSH URL is `git@github.com:brilliantlee2/flora_private.git`.
Patched dependencies are included in `vendor/`; no submodule setup is needed.

## Build environment

Run these commands from the repository root on Linux x86_64:

```bash
conda env create -f environment.yml
conda activate flora

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1
```

For an existing environment, use `conda env update -n flora -f environment.yml`
instead of creating it, then reactivate it. Use Python 3.11.x.
`environment.yml` includes the build tools and runtime dependencies;
`environment.runtime.yml` is for packaged executables, and `requirements.txt`
lists Python packages only. Both environment files include StringTie 2.2.2,
`gffread`, and `gffcompare`, required for transcript reconstruction and annotation.

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
`conda env config vars set LIBCLANG_PATH="$LIBCLANG_PATH"` and reactivate `flora`.

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

## Workflow output policy

Default light output retains the final tagged BAM/index, expression matrices,
cell/barcode mapping, RNA clustering, QC, `qc/metrics_summary.xlsx`, logs, and
the HTML report. After successful completion, reproducible read-level tables,
the pre-tagging aligned BAM, and assignment TSVs are removed. An empty
`alignment/` directory is also removed.

Additional outputs:

| Path | Contents |
| --- | --- |
| `transcriptome/<sample>.transcriptome.fa.gz` | Assembled transcript sequences. |
| `transcriptome/<sample>.transcriptome.gtf.gz` | Assembled transcript annotation. |
| `transcriptome/tool_versions.json` | Actual tool versions. |
| `transcriptome/assignment_summary.json` | Gene/transcript assignment status counts. |
| `matrix/<sample>.gene_raw_feature_bc_matrix/` | Gene-level integer UMI counts. |
| `matrix/<sample>.gene_processed_feature_bc_matrix/` | Filtered, normalized gene expression. |
| `matrix/<sample>.transcript_raw_feature_bc_matrix/` | Transcript-level integer UMI counts. |
| `matrix/<sample>.transcript_processed_feature_bc_matrix/` | Filtered, normalized transcript expression. |
| `matrix/<sample>.matrix_processing.json` | Matrix filtering settings and dimensions. |

Each MEX directory contains `matrix.mtx.gz`, `features.tsv.gz`, and
`barcodes.tsv.gz` (cell IDs). Processed matrices apply filters of at least 200
features/cell, 3 cells/feature and at most 20% mitochondrial counts, then normalize
to 10,000 counts/cell and apply `log10(1+x)`. Their cell set may differ from the
report's original cell calls. `--skip-isoform` omits transcript matrices.
The existing TSV matrices are retained.

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
