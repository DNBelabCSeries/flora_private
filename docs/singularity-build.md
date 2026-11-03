[English](singularity-build.md) | [中文](singularity-build_zh-CN.md)

# Build a Singularity Image

This guide was moved from the private README on 2026-10-27 to keep the main
build instructions short. It packages an existing Linux release; it does not
change Flora's analysis logic.

## Prerequisites

- Linux x86_64 and Singularity with working `--fakeroot` on the build host.
- A Flora release archive and the matching `packaging/singularity/Flora.def`.
- Network access for Ubuntu packages and Conda dependencies.
- Docker only if preparing the base image from a Docker image. An existing base
  SIF can be copied instead.

The examples match the current definition file, version 0.1.2. For another
release, update the definition's label and archive name, and the commands below.
The base used here is Ubuntu 22.04 (glibc 2.35). Check that the Flora binary's
`required_glibc` in `BUILD_INFO.txt` does not exceed the container's glibc.
The Flora Conda environment does not need to be active when building the SIF.

## 1. Prepare a build directory

From the private repository root, after building the release archive:

```bash
mkdir -p singularity_build
cp dist/Flora-0.1.2-linux-x86_64.tar.gz singularity_build/
cp packaging/singularity/Flora.def singularity_build/
cd singularity_build
```

## 2. Prepare the base image

If `flora-base-ubuntu22.04.sif` is already available, copy it into this directory
and skip this step. Otherwise, inspect and export the Docker image:

```bash
docker run --rm quay.io/nf-core/ubuntu:22.04 \
  bash -c 'cat /etc/os-release; getconf GNU_LIBC_VERSION'

docker save quay.io/nf-core/ubuntu:22.04 \
  -o flora-base-ubuntu22.04.tar

singularity build --fakeroot \
  flora-base-ubuntu22.04.sif \
  docker-archive://flora-base-ubuntu22.04.tar
```

The inspection should report Ubuntu 22.04 and glibc 2.35. If the Docker image is
not cached locally, Docker needs registry access. Do not make the Docker socket
world-writable to resolve permission errors; use an authorized Docker account
or ask the administrator. Importing an existing Docker archive does not require
access to the Docker daemon.

## 3. Download the local Miniforge installer

`Flora.def` copies this installer from the host; it does not download it during
the build. Reuse an existing trusted installer or download it once:

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

These four files must be in the build directory:

```text
Flora.def
Flora-0.1.2-linux-x86_64.tar.gz
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

## 4. Build and validate

```bash
env -u LD_LIBRARY_PATH \
singularity build --fakeroot \
  Flora-0.1.2-linux-x86_64.sif \
  Flora.def

singularity run --cleanenv Flora-0.1.2-linux-x86_64.sif --version
singularity run --cleanenv Flora-0.1.2-linux-x86_64.sif --help

sha256sum Flora-0.1.2-linux-x86_64.sif \
  > Flora-0.1.2-linux-x86_64.sif.sha256
sha256sum -c Flora-0.1.2-linux-x86_64.sif.sha256
```

Check that each command succeeds before continuing. The definition's `%test`
checks Python imports and required executables. Final image size depends on the
resolved dependencies. Cache directories are configured under `/tmp` by default.

If `--fakeroot` reports missing subuid mappings, ask the administrator to
configure the build host, or build on another suitable host and copy the final
SIF. If execution reports disabled user namespaces, the target host needs a
working administrator-configured Singularity runtime; changing Flora's files
does not fix that host restriction.

## 5. Run and distribute

Users need only the final SIF and a working Singularity installation, not Docker,
the base image, or a separate Flora Conda environment. Bind every data/reference/
output parent directory needed by the workflow. For example:

```bash
singularity run --cleanenv \
  --bind /data:/data \
  Flora-0.1.2-linux-x86_64.sif \
  --fastq-dir /data/sample_fastq_dir \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir /data/sample_output \
  --sample-id sample
```

Upload the final SIF and checksum as optional GitHub Release assets. Do not
commit images, Docker archives, or Miniforge installers to Git.

[Back to the private README](../README.md)
