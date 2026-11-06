#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash run_all.sh \
    --fastq /path/to/chunk_1.fq.gz [/path/to/chunk_2.fq.gz ...] \
    --tso-seq AAGACCGCTTGGCCTCCGACTTTCTGCG \
    --rtp-seq GAGGTCCATGAAGTGAGCATCTTCTGCG \
    --barcode-list-10bp /path/to/BC_1536.txt \
    --ref-dir /path/to/ref_dir \
    --out-dir /path/to/output_dir \
    [--gene-fasta /path/to/genome.fa] \
    [--junction-bed /path/to/genes.bed] \
    [--chrom-sizes /path/to/chrom_sizes.tsv] \
    [--gene-gtf /path/to/genes.gtf] \
    [--isoform-gtf /path/to/transcripts.gtf] \
    [--sample-id sample] \
    [--threads 32] \
    [--fastq-dir /path/to/fastq_chunks] \
    [--glycine-jobs 10] \
    [--glycine-threads 64] \
    [--glycine-outdir /path/to/glycine_out] \
    [--glycine-err 0.2,0.25] \
    [--glycine-shift 200,200] \
    [--min-len 300] \
    [--umi-len 41] \
    [--full-length-fastq /path/to/sample.full-length-plus-rescued.fq.gz] \
    [--skip-glycine] \
    [--skip-isoform] \
    [--upstream-only] \
    [--full-pipeline] \
    [--light-output] \
    [--full-output] \
    [--no-revcomp-whitelist] \
    [--cluster-threads 16] \
    [--exp-cells 5000] \
    [--min-q 2] \
    [--max-ed 2] \
    [--barcode-extract-mode fixed_seq|probe] \
    [--pair-min 10] \
    [--auto-pair-min-floor 10] \
    [--auto-pair-min-quantile 0.1] \
    [--top1-alpha 0.1] \
    [--top1-alpha-umi 0.3] \
    [--singlet-threshold 0.9] \
    [--dominance-min 0.8] \
    [--drop-umiA-ratio-gt 0.5] \
    [--gene-assign-mapq 30] \
    [--gene-assign-chunk-size 200000] \
    [--gene-tagging-backend rust|legacy|validate] \
    [--feature-backend ont|legacy] \
    [--ont-jobs 1] \
    [--ont-temp-dir /path/to/local_scratch] \
    [--stringtie-opts '-c 2'] \
    [--transcript-assign-mapq 30] \
    [--transcript-assign-chunk-size 200000] \
    [--ref-interval 1000] \
    [--cell-gene-max-reads 0] \
    [--save-merge-debug] \
    [--save-intermediate] \
    [--remove-final-bam] \
    [--require-pass-both-ends]

Notes:
  1. If --skip-glycine is not set, Glycine runs first and the pipeline uses:
     <glycine-outdir>/<sample-id>.full-length-plus-rescued.fq.gz
  2. Glycine is built into the Flora executable; no separate installation is needed.
     Multiple paths after --fastq run concurrently. --fastq-dir discovers only
     *.fastq.gz and *.fq.gz files in that directory (not subdirectories).
  3. --barcode-list-10bp accepts a plain text/CSV/TSV/XLSX file with 10bp barcodes
     and expands them into 26bp 3p/5p whitelists before Flora cell assignment.
  4. If --barcode-list-10bp is omitted, the bundled
     "BC_1536.txt" is used when present.
  5. --ref-dir looks for genome.fa, genes.bed, chrom_sizes.tsv, and genes.gtf.
  6. Explicit reference file arguments override files inferred from --ref-dir.
  7. --gene-fasta is the minimap2 reference FASTA. --genome-fa is accepted as an alias.
  8. If --isoform-gtf is omitted, --gene-gtf is reused for isoform assignment.
  9. If --pair-min is omitted, PAIR_MIN is auto-resolved as:
     max(--auto-pair-min-floor, quantile(pair support, --auto-pair-min-quantile)).
  10. Flora keeps the validated dual-end barcode model, but writes Sockeye-style
     BAM tags: CB/CR/UR plus dual-end custom tags C5/C3/U5/U3.
  11. --barcode-extract-mode fixed_seq is the current validated mode. probe is
     reserved for a future Sockeye-style local-alignment extractor and exits clearly.
  12. Light output is enabled by default. Trimmed cDNA reads stream directly to
      minimap2 without creating alignment/<sample>.cdna.fastq.gz on disk.
  13. Use --full-output or --save-intermediate to retain the trimmed cDNA FASTQ
      and other intermediate outputs for debugging.
     Use --remove-final-bam to remove every alignment/matrix BAM and index after
     matrices, QC, barnyard QC, and the report finish successfully.
  14. This mixed-species entrypoint expects a combined reference with prefixed
      contigs/IDs/names, for example hs_chr1 and mm_chr1.
  15. --skip-isoform skips isoform matrices (and legacy transcript assignment).
      Gene-level expression, barnyard QC, RNA QC, and the HTML report still run.
  16. Advanced knobs are available for assignment and memory control:
      --gene-assign-mapq, --gene-assign-chunk-size,
      --transcript-assign-mapq, --transcript-assign-chunk-size,
      --ref-interval, and --cell-gene-max-reads.
  17. --feature-backend ont (default) uses wf-single-cell v3.3.0-style sample
      transcriptome assembly and joint gene/transcript assignment. Requires
      stringtie, gffread and gffcompare. legacy uses the original GTF overlap.
      --gene-assign-mapq gates genes; --transcript-assign-mapq is legacy-only.
      ONT defaults to no per-cell/gene read cap (0); legacy defaults to 20000.
      ONT 3.3.0 was withdrawn upstream. This comparison backend pins historical
      source d46f92121c79c28a02105cf2444029e91dc4b609 and StringTie 2.2.2.
EOF
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] Command not found: $cmd" >&2
    exit 1
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing ${label}: $path" >&2
    exit 1
  fi
}

python_asset() {
  local source_path="$1"
  if [[ -f "${source_path}" ]]; then
    printf '%s\n' "${source_path}"
  elif [[ -f "${source_path}c" ]]; then
    printf '%s\n' "${source_path}c"
  else
    echo "[ERROR] Missing Python runtime asset: ${source_path} or ${source_path}c" >&2
    return 1
  fi
}

log() {
  echo "[run_all] $*"
}

SCRIPT_START_TS=$(date +%s)
STEP_TS=0

step_start() {
  STEP_TS=$(date +%s)
}

step_end() {
  local step_name="$1"
  local step_end_ts
  step_end_ts=$(date +%s)
  log "${step_name} elapsed: $((step_end_ts - STEP_TS))s"
}

write_parameters_tsv() {
  local out_tsv="$1"
  {
    printf "Parameter\tValue\n"
    printf "sample_id\t%s\n" "${SAMPLE_ID}"
    printf "skip_glycine\t%s\n" "${SKIP_GLYCINE}"
    printf "skip_isoform\t%s\n" "${SKIP_ISOFORM}"
    printf "upstream_only\t%s\n" "${UPSTREAM_ONLY}"
    printf "light_output\t%s\n" "${LIGHT_OUTPUT}"
    printf "remove_final_bam\t%s\n" "${REMOVE_FINAL_BAM}"
    printf "skip_matched_fastq\t%s\n" "${SKIP_MATCHED_FASTQ}"
    printf "skip_unmatched_fastq\t%s\n" "${SKIP_UNMATCHED_FASTQ}"
    printf "skip_cell_fastq\t%s\n" "${SKIP_CELL_FASTQ}"
    printf "fastq\t%s\n" "${FASTQ}"
    printf "full_length_fastq\t%s\n" "${FULL_LENGTH_FASTQ}"
    printf "alignment_input\t%s\n" "${ALIGNMENT_INPUT_MODE}"
    printf "glycine_outdir\t%s\n" "${GLYCINE_OUTDIR}"
    printf "barcode_list_10bp\t%s\n" "${BARCODE_LIST_10BP}"
    printf "tso_seq\t%s\n" "${TSO_SEQ}"
    printf "rtp_seq\t%s\n" "${RTP_SEQ}"
    printf "ref_dir\t%s\n" "${REF_DIR}"
    printf "gene_fasta\t%s\n" "${GENE_FASTA}"
    printf "junction_bed\t%s\n" "${JUNCTION_BED}"
    printf "chrom_sizes\t%s\n" "${CHROM_SIZES}"
    printf "gene_gtf\t%s\n" "${GENE_GTF}"
    printf "isoform_gtf\t%s\n" "${ISOFORM_GTF}"
    printf "out_dir\t%s\n" "${OUT_DIR}"
    printf "threads\t%s\n" "${THREADS}"
    printf "cluster_threads\t%s\n" "${CLUSTER_THREADS}"
    printf "glycine_jobs\t%s\n" "${GLYCINE_JOBS}"
    printf "glycine_threads\t%s\n" "${GLYCINE_THREADS}"
    printf "exp_cells\t%s\n" "${EXP_CELLS}"
    printf "min_q\t%s\n" "${MIN_Q}"
    printf "max_ed\t%s\n" "${MAX_ED}"
    printf "barcode_extract_mode\t%s\n" "${BARCODE_EXTRACT_MODE}"
    printf "pair_min\t%s\n" "${PAIR_MIN}"
    printf "auto_pair_min_floor\t%s\n" "${AUTO_PAIR_MIN_FLOOR}"
    printf "auto_pair_min_quantile\t%s\n" "${AUTO_PAIR_MIN_QUANTILE}"
    printf "top1_alpha\t%s\n" "${TOP1_ALPHA}"
    printf "top1_alpha_umi\t%s\n" "${TOP1_ALPHA_UMI}"
    printf "singlet_threshold\t%s\n" "${SINGLET_THRESHOLD}"
    printf "dominance_min\t%s\n" "${DOMINANCE_MIN}"
    printf "drop_umiA_ratio_gt\t%s\n" "${DROP_UMIA_RATIO_GT}"
    printf "require_pass_both_ends\t%s\n" "${REQUIRE_PASS_BOTH_ENDS}"
    printf "glycine_err\t%s\n" "${GLYCINE_ERR}"
    printf "glycine_shift\t%s\n" "${GLYCINE_SHIFT}"
    printf "min_len\t%s\n" "${MIN_LEN}"
    printf "umi_len\t%s\n" "${UMI_LEN}"
    printf "gene_assign_mapq\t%s\n" "${GENE_ASSIGN_MAPQ}"
    printf "gene_assign_chunk_size\t%s\n" "${GENE_ASSIGN_CHUNK_SIZE}"
    printf "gene_tagging_backend\t%s\n" "${GENE_TAGGING_BACKEND}"
    printf "feature_backend\t%s\n" "${FEATURE_BACKEND}"
    if [[ "${FEATURE_BACKEND}" == "ont" ]]; then
      printf "feature_reference\twf-single-cell v3.3.0\n"
      printf "feature_reference_commit\td46f92121c79c28a02105cf2444029e91dc4b609\n"
      printf "feature_reference_withdrawn\ttrue\n"
    else
      printf "feature_reference\tFlora GTF-overlap\n"
    fi
    printf "stringtie_opts\t%s\n" "${STRINGTIE_OPTS}"
    printf "ont_jobs\t%s\n" "${ONT_JOBS}"
    printf "ont_temp_dir\t%s\n" "${ONT_TEMP_DIR:-system temp/TMPDIR}"
    printf "transcript_assign_mapq\t%s\n" "${TRANSCRIPT_ASSIGN_MAPQ}"
    printf "transcript_assign_chunk_size\t%s\n" "${TRANSCRIPT_ASSIGN_CHUNK_SIZE}"
    printf "ref_interval\t%s\n" "${REF_INTERVAL}"
    printf "cell_gene_max_reads\t%s\n" "${CELL_GENE_MAX_READS}"
  } > "${out_tsv}"
}

resolve_exec() {
  local rust_name="$1"
  local py_path="$2"
  local rust_candidates=(
    "${SCRIPT_DIR}/target/release/${rust_name}"
  )
  local rust_bin=""
  for candidate in "${rust_candidates[@]}"; do
    if [[ -x "${candidate}" ]]; then
      rust_bin="${candidate}"
      break
    fi
  done
  if [[ -n "${rust_bin}" ]]; then
    echo "${rust_bin}"
    return 0
  fi
  py_path="$(python_asset "${py_path}")"
  echo "python3 ${py_path}"
}

run_stage() {
  local rust_name="$1"
  local py_path="$2"
  shift 2
  local rust_candidates=(
    "${SCRIPT_DIR}/target/release/${rust_name}"
  )
  local rust_bin=""
  local stage_start_ts stage_end_ts stage_status stage_impl
  stage_start_ts=$(date +%s)
  for candidate in "${rust_candidates[@]}"; do
    if [[ -x "${candidate}" ]]; then
      rust_bin="${candidate}"
      break
    fi
  done
  if [[ -n "${rust_bin}" ]]; then
    stage_impl="rust"
    if "${rust_bin}" "$@"; then stage_status=0; else stage_status=$?; fi
  else
    stage_impl="python"
    py_path="$(python_asset "${py_path}")"
    if python3 "${py_path}" "$@"; then stage_status=0; else stage_status=$?; fi
  fi
  stage_end_ts=$(date +%s)
  log "Stage ${rust_name} (${stage_impl}) elapsed: $((stage_end_ts - stage_start_ts))s"
  return "${stage_status}"
}

abspath_for_output() {
  local path="$1"
  mkdir -p "$(dirname "$path")"
  echo "$(cd "$(dirname "$path")" && pwd)/$(basename "$path")"
}

remove_intermediate() {
  local path="$1"
  if [[ -e "${path}" || -L "${path}" ]]; then
    rm -f "${path}"
    log "Removed intermediate: ${path}"
  fi
}

run_alignment_pipeline() {
  local input_fastq="$1"
  minimap2 -ax splice -uf --MD -I 16G -t "${THREADS}" -K 10M \
    --junc-bed "${JUNCTION_BED}" \
    --cap-kalloc 100m \
    "${GENE_FASTA}" "${input_fastq}" \
  | samtools view -uh --no-PG - \
  | samtools view -uh -F 256 - \
  | samtools sort --no-PG -@ "${THREADS}" -o "${ALIGNED_BAM}" -
}

cleanup_alignment_stream() {
  if [[ "${ALIGNMENT_STREAM_ACTIVE:-0}" -ne 1 ]]; then
    return
  fi
  { exec 9>&-; } 2>/dev/null || true
  if [[ -n "${ALIGNMENT_STREAM_PID:-}" ]]; then
    kill "${ALIGNMENT_STREAM_PID}" 2>/dev/null || true
    wait "${ALIGNMENT_STREAM_PID}" 2>/dev/null || true
  fi
  if [[ -n "${ALIGNMENT_FASTQ:-}" && -p "${ALIGNMENT_FASTQ}" ]]; then
    rm -f "${ALIGNMENT_FASTQ}"
  fi
}

cleanup_large_intermediates() {
  local path
  local paths=(
    "${ALIGN_DIR}/${SAMPLE_ID}.cdna.fastq.gz"
    "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
    "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam.bai"
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam"
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam.bai"
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.bed"
    "${ALIGN_DIR}/${SAMPLE_ID}.read_tags.tsv"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam.bai"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.read_gene_assigns.tsv"
    "${MATRIX_DIR}/${SAMPLE_ID}.read_transcript_assigns.tsv"
    "${MATRIX_DIR}/${SAMPLE_ID}.cell_umi_gene.tsv"
    "${UPSTREAM_DIR}/ReadIDs_UMI_BC_clean.csv"
    "${UPSTREAM_DIR}/read_assigned_cell.csv"
    "${UPSTREAM_DIR}/correction_map_3p.tsv"
    "${UPSTREAM_DIR}/correction_map_5p.tsv"
    "${UPSTREAM_DIR}/read_qc_summary.json"
    "${UPSTREAM_DIR}/full_length_fastq_count.txt"
    "${QC_DIR}/cell_umi_gene.tsv"
    "${QC_DIR}/filtered.sorted.bam"
    "${QC_DIR}/filtered.sorted.bam.bai"
  )
  for path in "${paths[@]}"; do
    remove_intermediate "${path}"
  done
  rmdir "${ALIGN_DIR}" 2>/dev/null || true
}

cleanup_all_bam_outputs() {
  local path
  local paths=(
    "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
    "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam.bai"
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam"
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam.bai"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam.bai"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.tagged.sorted.bam"
    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.tagged.sorted.bam.bai"
    "${QC_DIR}/filtered.sorted.bam"
    "${QC_DIR}/filtered.sorted.bam.bai"
  )
  for path in "${paths[@]}"; do
    remove_intermediate "${path}"
  done
  rmdir "${ALIGN_DIR}" 2>/dev/null || true
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNSTREAM_DIR="${SCRIPT_DIR}/scripts"

FASTQ=""
FASTQ_INPUTS=()
FASTQ_DIR=""
FULL_LENGTH_FASTQ=""
TSO_SEQ="AAGACCGCTTGGCCTCCGACTTTCTGCG"
RTP_SEQ="GAGGTCCATGAAGTGAGCATCTTCTGCG"
BARCODE_LIST_10BP=""
REF_DIR=""
GENE_FASTA=""
JUNCTION_BED=""
CHROM_SIZES=""
GENE_GTF=""
ISOFORM_GTF=""
OUT_DIR=""
SAMPLE_ID="sample"
THREADS=32
CLUSTER_THREADS=16
GLYCINE_JOBS=10
GLYCINE_THREADS=64
EXP_CELLS=5000
MIN_Q=2
MAX_ED=2
BARCODE_EXTRACT_MODE="fixed_seq"
PAIR_MIN=""
AUTO_PAIR_MIN_FLOOR=10
AUTO_PAIR_MIN_QUANTILE=0.1
TOP1_ALPHA=0.1
TOP1_ALPHA_UMI=0.3
SINGLET_THRESHOLD=0.9
DOMINANCE_MIN=0.8
DROP_UMIA_RATIO_GT=0.5
GLYCINE_OUTDIR=""
GLYCINE_ERR="0.2,0.25"
GLYCINE_SHIFT="200,200"
MIN_LEN=300
UMI_LEN=41
SKIP_GLYCINE=0
SKIP_ISOFORM=0
UPSTREAM_ONLY=0
LIGHT_OUTPUT=1
SKIP_MATCHED_FASTQ=0
SKIP_UNMATCHED_FASTQ=0
SKIP_CELL_FASTQ=0
REVCOMP_WHITELIST=1
SAVE_MERGE_DEBUG=0
SAVE_INTERMEDIATE=0
REMOVE_FINAL_BAM=0
REQUIRE_PASS_BOTH_ENDS=0
GENE_ASSIGN_MAPQ=30
GENE_ASSIGN_CHUNK_SIZE=200000
GENE_TAGGING_BACKEND="rust"
FEATURE_BACKEND="ont"
ONT_JOBS=1
ONT_TEMP_DIR=""
STRINGTIE_OPTS="-c 2"
TRANSCRIPT_ASSIGN_MAPQ=30
TRANSCRIPT_ASSIGN_CHUNK_SIZE=200000
REF_INTERVAL=1000
CELL_GENE_MAX_READS=20000
CELL_GENE_MAX_READS_EXPLICIT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        FASTQ_INPUTS+=("$1")
        shift
      done
      ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --full-length-fastq) FULL_LENGTH_FASTQ="$2"; shift 2 ;;
    --tso-seq|--tso_seq) TSO_SEQ="$2"; shift 2 ;;
    --rtp-seq|--rtp_seq) RTP_SEQ="$2"; shift 2 ;;
    --barcode-list-10bp|--barcode_list_10bp) BARCODE_LIST_10BP="$2"; shift 2 ;;
    --whitelist-3p|--whitelist-5p)
      echo "[ERROR] Flora does not accept --whitelist-3p/--whitelist-5p. Use --barcode-list-10bp instead." >&2
      exit 1
      ;;
    --ref-dir|--ref_dir) REF_DIR="$2"; shift 2 ;;
    --gene-fasta|--genome-fa) GENE_FASTA="$2"; shift 2 ;;
    --junction-bed) JUNCTION_BED="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --gene-gtf) GENE_GTF="$2"; shift 2 ;;
    --isoform-gtf) ISOFORM_GTF="$2"; shift 2 ;;
    --out-dir|--outdir) OUT_DIR="$2"; shift 2 ;;
    --sample-id|--sample) SAMPLE_ID="$2"; shift 2 ;;
    --threads|--thread) THREADS="$2"; shift 2 ;;
    --cluster-threads) CLUSTER_THREADS="$2"; shift 2 ;;
    --glycine-jobs) GLYCINE_JOBS="$2"; shift 2 ;;
    --glycine-threads) GLYCINE_THREADS="$2"; shift 2 ;;
    --exp-cells) EXP_CELLS="$2"; shift 2 ;;
    --min-q) MIN_Q="$2"; shift 2 ;;
    --max-ed) MAX_ED="$2"; shift 2 ;;
    --barcode-extract-mode|--barcode_extract_mode) BARCODE_EXTRACT_MODE="$2"; shift 2 ;;
    --pair-min) PAIR_MIN="$2"; shift 2 ;;
    --auto-pair-min-floor) AUTO_PAIR_MIN_FLOOR="$2"; shift 2 ;;
    --auto-pair-min-quantile) AUTO_PAIR_MIN_QUANTILE="$2"; shift 2 ;;
    --top1-alpha) TOP1_ALPHA="$2"; shift 2 ;;
    --top1-alpha-umi) TOP1_ALPHA_UMI="$2"; shift 2 ;;
    --singlet-threshold) SINGLET_THRESHOLD="$2"; shift 2 ;;
    --dominance-min) DOMINANCE_MIN="$2"; shift 2 ;;
    --drop-umiA-ratio-gt) DROP_UMIA_RATIO_GT="$2"; shift 2 ;;
    --glycine-outdir) GLYCINE_OUTDIR="$2"; shift 2 ;;
    --glycine-err|--err) GLYCINE_ERR="$2"; shift 2 ;;
    --glycine-shift|--shift) GLYCINE_SHIFT="$2"; shift 2 ;;
    --min-len|--min_len) MIN_LEN="$2"; shift 2 ;;
    --umi-len|--umi_len) UMI_LEN="$2"; shift 2 ;;
    --skip-glycine) SKIP_GLYCINE=1; shift ;;
    --skip-isoform) SKIP_ISOFORM=1; shift ;;
    --upstream-only) UPSTREAM_ONLY=1; shift ;;
    --full-pipeline) UPSTREAM_ONLY=0; shift ;;
    --light-output) LIGHT_OUTPUT=1; shift ;;
    --full-output) LIGHT_OUTPUT=0; shift ;;
    --no-revcomp-whitelist) REVCOMP_WHITELIST=0; shift ;;
    --save-merge-debug) SAVE_MERGE_DEBUG=1; shift ;;
    --save-intermediate) SAVE_INTERMEDIATE=1; shift ;;
    --remove-final-bam) REMOVE_FINAL_BAM=1; shift ;;
    --require-pass-both-ends) REQUIRE_PASS_BOTH_ENDS=1; shift ;;
    --gene-assign-mapq) GENE_ASSIGN_MAPQ="$2"; shift 2 ;;
    --gene-assign-chunk-size) GENE_ASSIGN_CHUNK_SIZE="$2"; shift 2 ;;
    --gene-tagging-backend) GENE_TAGGING_BACKEND="$2"; shift 2 ;;
    --feature-backend) FEATURE_BACKEND="$2"; shift 2 ;;
    --ont-jobs) ONT_JOBS="$2"; shift 2 ;;
    --ont-temp-dir) ONT_TEMP_DIR="$2"; shift 2 ;;
    --stringtie-opts) STRINGTIE_OPTS="$2"; shift 2 ;;
    --stringtie-opts=*) STRINGTIE_OPTS="${1#*=}"; shift ;;
    --transcript-assign-mapq) TRANSCRIPT_ASSIGN_MAPQ="$2"; shift 2 ;;
    --transcript-assign-chunk-size) TRANSCRIPT_ASSIGN_CHUNK_SIZE="$2"; shift 2 ;;
    --ref-interval) REF_INTERVAL="$2"; shift 2 ;;
    --cell-gene-max-reads) CELL_GENE_MAX_READS="$2"; CELL_GENE_MAX_READS_EXPLICIT=1; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ! "${ONT_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] --ont-jobs must be a positive integer." >&2
  exit 1
fi
if [[ "${FEATURE_BACKEND}" != "ont" && "${FEATURE_BACKEND}" != "legacy" ]]; then
  echo "[ERROR] --feature-backend must be ont or legacy." >&2
  exit 1
fi
if [[ "${FEATURE_BACKEND}" == "ont" && "${CELL_GENE_MAX_READS_EXPLICIT}" -eq 0 ]]; then
  CELL_GENE_MAX_READS=0
fi

if [[ "${GENE_TAGGING_BACKEND}" != "rust" && "${GENE_TAGGING_BACKEND}" != "legacy" && "${GENE_TAGGING_BACKEND}" != "validate" ]]; then
  echo "[ERROR] --gene-tagging-backend must be rust, legacy, or validate." >&2
  exit 1
fi

if [[ "${LIGHT_OUTPUT}" -eq 1 ]]; then
  SKIP_MATCHED_FASTQ=1
  SKIP_UNMATCHED_FASTQ=1
  SKIP_CELL_FASTQ=1
else
  SKIP_MATCHED_FASTQ=0
  SKIP_UNMATCHED_FASTQ=0
  SKIP_CELL_FASTQ=0
fi

if [[ "${BARCODE_EXTRACT_MODE}" != "fixed_seq" && "${BARCODE_EXTRACT_MODE}" != "probe" ]]; then
  echo "[ERROR] --barcode-extract-mode must be fixed_seq or probe, got: ${BARCODE_EXTRACT_MODE}" >&2
  exit 1
fi

if [[ -n "${REF_DIR}" ]]; then
  if [[ ! -d "${REF_DIR}" ]]; then
    echo "[ERROR] Missing reference directory: ${REF_DIR}" >&2
    exit 1
  fi
  REF_DIR="$(cd "${REF_DIR}" && pwd)"
  [[ -z "${GENE_FASTA}" ]] && GENE_FASTA="${REF_DIR}/genome.fa"
  [[ -z "${JUNCTION_BED}" ]] && JUNCTION_BED="${REF_DIR}/genes.bed"
  [[ -z "${CHROM_SIZES}" ]] && CHROM_SIZES="${REF_DIR}/chrom_sizes.tsv"
  [[ -z "${GENE_GTF}" ]] && GENE_GTF="${REF_DIR}/genes.gtf"
fi

if [[ -z "${OUT_DIR}" || -z "${GENE_FASTA}" || -z "${JUNCTION_BED}" || -z "${CHROM_SIZES}" || -z "${GENE_GTF}" ]]; then
  echo "[ERROR] Missing required arguments." >&2
  usage
  exit 1
fi

if [[ -z "${ISOFORM_GTF}" ]]; then
  ISOFORM_GTF="${GENE_GTF}"
fi

OUT_DIR="$(abspath_for_output "${OUT_DIR}")"
if [[ -n "${ONT_TEMP_DIR}" ]]; then
  ONT_TEMP_DIR="$(abspath_for_output "${ONT_TEMP_DIR}")"
fi
BARCODE_DIR="${OUT_DIR}/barcodes"
require_cmd python3
if [[ -z "${BARCODE_LIST_10BP}" ]]; then
  DEFAULT_BARCODE_LIST="${SCRIPT_DIR}/BC_1536.txt"
  if [[ -f "${DEFAULT_BARCODE_LIST}" ]]; then
    BARCODE_LIST_10BP="${DEFAULT_BARCODE_LIST}"
  fi
fi
if [[ -z "${BARCODE_LIST_10BP}" ]]; then
  echo "[ERROR] Missing --barcode-list-10bp and bundled BC_1536.txt was not found." >&2
  usage
  exit 1
fi
mkdir -p "${BARCODE_DIR}"
WHITELIST_3P="${BARCODE_DIR}/barcode_3prime_26bp.txt"
WHITELIST_5P="${BARCODE_DIR}/barcode_5prime_26bp.txt"
log "Generating 26bp barcode whitelists from 10bp list: ${BARCODE_LIST_10BP}"
require_file "${BARCODE_LIST_10BP}" "10bp barcode list"
run_stage generate_26bp_whitelists "${DOWNSTREAM_DIR}/generate_26bp_whitelists.py" \
  --barcode-list-10bp "${BARCODE_LIST_10BP}" \
  --out-3p "${WHITELIST_3P}" \
  --out-5p "${WHITELIST_5P}" 2>&1 | tee "${OUT_DIR}/barcode_whitelist_generation.log"
if [[ "${SKIP_GLYCINE}" -eq 0 ]]; then
  if [[ -z "${GLYCINE_OUTDIR}" ]]; then
    GLYCINE_OUTDIR="${OUT_DIR}/glycine"
  else
    GLYCINE_OUTDIR="$(abspath_for_output "${GLYCINE_OUTDIR}")"
  fi
else
  GLYCINE_OUTDIR=""
fi

UPSTREAM_DIR="${OUT_DIR}/upstream"
ALIGN_DIR="${OUT_DIR}/alignment"
MATRIX_DIR="${OUT_DIR}/matrix"
QC_DIR="${OUT_DIR}/qc"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${UPSTREAM_DIR}" "${ALIGN_DIR}" "${MATRIX_DIR}" "${QC_DIR}" "${LOG_DIR}"
if [[ "${SKIP_GLYCINE}" -eq 0 ]]; then
  mkdir -p "${GLYCINE_OUTDIR}"
fi

require_cmd python3
require_cmd samtools
require_cmd minimap2
if [[ "${FEATURE_BACKEND}" == "legacy" && ( "${GENE_TAGGING_BACKEND}" != "rust" || "${SKIP_ISOFORM}" -eq 0 ) ]]; then
  require_cmd bedtools
fi
if [[ "${FEATURE_BACKEND}" == "ont" && "${UPSTREAM_ONLY}" -eq 0 ]]; then
  require_cmd stringtie
  require_cmd gffread
  require_cmd gffcompare
  if [[ "$(stringtie --version)" != "2.2.2" ]]; then
    echo "[ERROR] ONT 3.3.0 comparison requires StringTie 2.2.2; update environment.yml." >&2
    exit 1
  fi
  python_asset "${DOWNSTREAM_DIR}/ont_features.py" >/dev/null
  python_asset "${DOWNSTREAM_DIR}/ont_matrices.py" >/dev/null
fi
require_file "${WHITELIST_3P}" "3p whitelist"
require_file "${WHITELIST_5P}" "5p whitelist"
require_file "${GENE_FASTA}" "gene FASTA"
require_file "${JUNCTION_BED}" "junction BED"
require_file "${CHROM_SIZES}" "chrom sizes"
require_file "${GENE_GTF}" "gene GTF"
require_file "${ISOFORM_GTF}" "isoform GTF"
python_asset "${SCRIPT_DIR}/main.py" >/dev/null

for script in prepare_read_tags.py add_cb_ur_tags.py assign_genes.py add_gene_tags.py cluster_umis_allbam.py cell_umi_gene_table.py gene_expression.py rna_cluster_analysis.py assign_transcripts.py isoform_expression.py rna_qc_metrics_mixed.py Saturation.py read_qc_summary.py metrics_summary.py build_report.py barnyard_qc.py generate_26bp_whitelists.py build_mixed_species_reference.py generate_knee_plots.py; do
  python_asset "${DOWNSTREAM_DIR}/${script}" >/dev/null
done
require_file "${DOWNSTREAM_DIR}/report_template.html" "report_template.html"
require_file "${DOWNSTREAM_DIR}/plotly-2.26.0.min.js" "plotly-2.26.0.min.js"

if [[ -f "${SCRIPT_DIR}/Cargo.toml" ]] && command -v cargo >/dev/null 2>&1; then
  log "Flora Rust project detected, attempting cargo build --release"
  (
    cd "${SCRIPT_DIR}"
    cargo build --release
  ) || log "cargo build failed; Python fallbacks remain available for supported downstream stages"
fi

FLORA_BIN="${SCRIPT_DIR}/target/release/flora"
if [[ "${SKIP_GLYCINE}" -eq 0 && ! -x "${FLORA_BIN}" ]]; then
  echo "[ERROR] Flora executable not found: ${FLORA_BIN}" >&2
  echo "[ERROR] Run 'cargo build --release' in ${SCRIPT_DIR} before using integrated Glycine." >&2
  exit 1
fi

if [[ "${SKIP_GLYCINE}" -eq 0 ]]; then
  if [[ ${#FASTQ_INPUTS[@]} -gt 0 && -n "${FASTQ_DIR}" ]]; then
    echo "[ERROR] --fastq and --fastq-dir cannot be used together." >&2
    exit 1
  fi
  if [[ ${#FASTQ_INPUTS[@]} -eq 0 && -z "${FASTQ_DIR}" ]] || [[ -z "${TSO_SEQ}" || -z "${RTP_SEQ}" ]]; then
    echo "[ERROR] --fastq, --tso-seq, and --rtp-seq are required unless --skip-glycine is set." >&2
    exit 1
  fi
  for fastq_path in "${FASTQ_INPUTS[@]}"; do
    require_file "${fastq_path}" "raw FASTQ"
  done
  if [[ -n "${FASTQ_DIR}" && ! -d "${FASTQ_DIR}" ]]; then
    echo "[ERROR] FASTQ directory does not exist: ${FASTQ_DIR}" >&2
    exit 1
  fi
  if [[ ${#FASTQ_INPUTS[@]} -gt 0 ]]; then
    FASTQ="$(IFS=';'; echo "${FASTQ_INPUTS[*]}")"
    GLYCINE_FASTQ_ARGS=(--fastq "${FASTQ_INPUTS[@]}")
  else
    FASTQ="${FASTQ_DIR}"
    GLYCINE_FASTQ_ARGS=(--fastq-dir "${FASTQ_DIR}")
  fi
  step_start
  log "Step 0/6: running parallel Glycine full-length identification"
  GLYCINE_STAGE_TS=$(date +%s)
  if "${FLORA_BIN}" glycine \
    "${GLYCINE_FASTQ_ARGS[@]}" \
    --tso_seq "${TSO_SEQ}" \
    --rtp_seq "${RTP_SEQ}" \
    --outdir "${GLYCINE_OUTDIR}" \
    --err "${GLYCINE_ERR}" \
    --min_len "${MIN_LEN}" \
    --shift "${GLYCINE_SHIFT}" \
    --umi_len "${UMI_LEN}" \
    --sample "${SAMPLE_ID}" \
    --jobs "${GLYCINE_JOBS}" \
    --total-threads "${GLYCINE_THREADS}" 2>&1 | tee "${LOG_DIR}/00_glycine.log"; then
    GLYCINE_STAGE_STATUS=0
  else
    GLYCINE_STAGE_STATUS=$?
  fi
  log "Stage glycine (rust) elapsed: $(($(date +%s) - GLYCINE_STAGE_TS))s"
  if (( GLYCINE_STAGE_STATUS != 0 )); then exit "${GLYCINE_STAGE_STATUS}"; fi
  step_end "Step 0/6"

  FULL_LENGTH_FASTQ="${GLYCINE_OUTDIR}/${SAMPLE_ID}.full-length-plus-rescued.fq.gz"
  GLYCINE_STATS="${GLYCINE_OUTDIR}/${SAMPLE_ID}.identifying_statistic.txt"
  require_file "${GLYCINE_STATS}" "Glycine identifying statistics"
else
  if [[ -z "${FULL_LENGTH_FASTQ}" ]]; then
    echo "[ERROR] --full-length-fastq is required with --skip-glycine." >&2
    exit 1
  fi
fi
require_file "${FULL_LENGTH_FASTQ}" "full-length FASTQ"

ALIGN_LOG="${LOG_DIR}/01_alignment.log"
UPSTREAM_LOG="${LOG_DIR}/02_dual_barcode.log"
TAG_LOG="${LOG_DIR}/03_bam_tagging.log"
GENE_LOG="${LOG_DIR}/04_gene.log"
ISOFORM_LOG="${LOG_DIR}/05_isoform.log"
QC_LOG="${LOG_DIR}/06_qc.log"
ALIGNED_BAM="${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
ALIGNMENT_STREAMING=0
ALIGNMENT_STREAM_ACTIVE=0
ALIGNMENT_STREAM_PID=""
ALIGNMENT_STAGE_TS=0
ALIGNMENT_INPUT_MODE="persisted barcode/UMI-trimmed cDNA FASTQ"
if [[ "${LIGHT_OUTPUT}" -eq 1 && "${SAVE_INTERMEDIATE}" -eq 0 ]]; then
  ALIGNMENT_STREAMING=1
  ALIGNMENT_INPUT_MODE="streamed barcode/UMI-trimmed cDNA FASTQ"
  ALIGNMENT_FASTQ="${ALIGN_DIR}/${SAMPLE_ID}.cdna.fastq.pipe"
else
  ALIGNMENT_FASTQ="${ALIGN_DIR}/${SAMPLE_ID}.cdna.fastq.gz"
fi

if [[ "${UPSTREAM_ONLY}" -eq 0 && "${ALIGNMENT_STREAMING}" -eq 1 ]]; then
  rm -f "${ALIGNMENT_FASTQ}"
  mkfifo "${ALIGNMENT_FASTQ}"
  exec 9<>"${ALIGNMENT_FASTQ}"
  ALIGNMENT_STREAM_ACTIVE=1
  trap cleanup_alignment_stream EXIT
  ALIGNMENT_STAGE_TS=$(date +%s)
  log "Streaming trimmed cDNA reads directly to minimap2 (no alignment FASTQ on disk)"
  (
    exec 9>&-
    run_alignment_pipeline "${ALIGNMENT_FASTQ}" 2>&1 | tee "${ALIGN_LOG}"
  ) &
  ALIGNMENT_STREAM_PID=$!
fi

step_start
log "Step 1/6: Flora barcode/UMI extraction, cDNA trimming, and cell assignment"
MAIN_ARGS=(
  "${FULL_LENGTH_FASTQ}"
  --full-bc-whitelist-3p "${WHITELIST_3P}"
  --full-bc-whitelist-5p "${WHITELIST_5P}"
  --out_dir "${UPSTREAM_DIR}"
  --threads "${THREADS}"
  --batch_size 100000
  --assign_batchsize 10000
  --exp_cells "${EXP_CELLS}"
  --minQ "${MIN_Q}"
  --max_ed "${MAX_ED}"
  --barcode_extract_mode "${BARCODE_EXTRACT_MODE}"
  --auto_pair_min_floor "${AUTO_PAIR_MIN_FLOOR}"
  --auto_pair_min_quantile "${AUTO_PAIR_MIN_QUANTILE}"
  --TOP1_ALPHA "${TOP1_ALPHA}"
  --TOP1_ALPHA_UMI "${TOP1_ALPHA_UMI}"
  --dominance_min "${DOMINANCE_MIN}"
  --drop_umiA_ratio_gt "${DROP_UMIA_RATIO_GT}"
)
if [[ "${UPSTREAM_ONLY}" -eq 0 ]]; then
  MAIN_ARGS+=(--alignment-fastq "${ALIGNMENT_FASTQ}")
fi
if [[ -n "${PAIR_MIN}" ]]; then
  MAIN_ARGS+=(--PAIR_MIN "${PAIR_MIN}")
fi
if [[ "${REVCOMP_WHITELIST}" -eq 0 ]]; then
  MAIN_ARGS+=(--no-revcomp-whitelist)
fi
if [[ "${SAVE_MERGE_DEBUG}" -eq 1 ]]; then
  MAIN_ARGS+=(--save_merge_debug)
fi
if [[ "${SAVE_INTERMEDIATE}" -eq 1 ]]; then
  MAIN_ARGS+=(--save-intermediate)
fi
if [[ "${LIGHT_OUTPUT}" -eq 1 ]]; then
  MAIN_ARGS+=(--light-output)
else
  MAIN_ARGS+=(--full-output)
fi
if [[ "${REQUIRE_PASS_BOTH_ENDS}" -eq 1 ]]; then
  MAIN_ARGS+=(--require_pass_both_ends)
fi
FLORA_TRIM_CACHE_DIR="${ONT_TEMP_DIR}" run_stage flora "${SCRIPT_DIR}/main.py" "${MAIN_ARGS[@]}" 2>&1 | tee "${UPSTREAM_LOG}"
if [[ "${ALIGNMENT_STREAM_ACTIVE}" -eq 1 ]]; then
  exec 9>&-
fi
KNEE_ARGS=(
  --counts-3p "${UPSTREAM_DIR}/barcode_counts_3p.tsv"
  --counts-5p "${UPSTREAM_DIR}/barcode_counts_5p.tsv"
  --full-whitelist-3p "${WHITELIST_3P}"
  --full-whitelist-5p "${WHITELIST_5P}"
  --exp-cells "${EXP_CELLS}"
  --out-3p "${UPSTREAM_DIR}/knee_plot_3p.png"
  --out-5p "${UPSTREAM_DIR}/knee_plot_5p.png"
)
if [[ "${REVCOMP_WHITELIST}" -eq 0 ]]; then
  KNEE_ARGS+=(--no-revcomp-whitelist)
fi
python3 "$(python_asset "${DOWNSTREAM_DIR}/generate_knee_plots.py")" "${KNEE_ARGS[@]}" 2>&1 | tee -a "${UPSTREAM_LOG}"
step_end "Step 1/6"

READ_ASSIGNED_CELL="${UPSTREAM_DIR}/read_assigned_cell.csv"
BARCODE_VALIDITY_SUMMARY="${UPSTREAM_DIR}/barcode_validity_summary.tsv"
CELL_READS_FASTQ="${UPSTREAM_DIR}/cell_reads.fastq.gz"
MATCHED_READS_FASTQ="${UPSTREAM_DIR}/matched_reads.fastq.gz"
require_file "${READ_ASSIGNED_CELL}" "read_assigned_cell.csv"
require_file "${BARCODE_VALIDITY_SUMMARY}" "barcode_validity_summary.tsv"
if [[ "${SKIP_CELL_FASTQ}" -eq 0 ]]; then
  require_file "${CELL_READS_FASTQ}" "cell_reads.fastq.gz"
fi
if [[ "${SKIP_MATCHED_FASTQ}" -eq 0 ]]; then
  require_file "${MATCHED_READS_FASTQ}" "matched_reads.fastq.gz"
fi

if [[ "${UPSTREAM_ONLY}" -eq 1 ]]; then
  log "Upstream-only mode enabled; stopping after barcode merge / cell assignment."
  echo
  echo "[run_all] Upstream outputs:"
  echo "  upstream dir              : ${UPSTREAM_DIR}"
  echo "  barcode validity summary  : ${BARCODE_VALIDITY_SUMMARY}"
  echo "  assign stats              : ${UPSTREAM_DIR}/assign_stats.tsv"
  echo "  assigned reads            : ${READ_ASSIGNED_CELL}"
  echo "  cell read stats           : ${UPSTREAM_DIR}/cell_read_stats.csv"
  echo "  pair counts kept          : ${UPSTREAM_DIR}/pair_counts_kept.csv"
  if [[ "${SAVE_MERGE_DEBUG}" -eq 1 ]]; then
    echo "  pair counts all           : ${UPSTREAM_DIR}/pair_counts_all.csv"
    echo "  pair counts pairmin kept  : ${UPSTREAM_DIR}/pair_counts_pairmin_kept.csv"
    echo "  dropped pairs             : ${UPSTREAM_DIR}/dropped_pairs.csv"
    echo "  core cells debug          : ${UPSTREAM_DIR}/core_cells_debug.csv"
    echo "  assigned reads debug      : ${UPSTREAM_DIR}/read_assigned_debug.csv"
  fi
  exit 0
fi

step_start
log "Step 2/6: ONT-style minimap2 alignment of barcode/UMI-trimmed cDNA reads"
if [[ "${ALIGNMENT_STREAMING}" -eq 1 ]]; then
  if wait "${ALIGNMENT_STREAM_PID}"; then
    ALIGNMENT_STAGE_STATUS=0
  else
    ALIGNMENT_STAGE_STATUS=$?
  fi
  ALIGNMENT_STREAM_PID=""
  rm -f "${ALIGNMENT_FASTQ}"
  ALIGNMENT_STREAM_ACTIVE=0
  trap - EXIT
else
  require_file "${ALIGNMENT_FASTQ}" "barcode/UMI-trimmed alignment FASTQ"
  ALIGNMENT_STAGE_TS=$(date +%s)
  if run_alignment_pipeline "${ALIGNMENT_FASTQ}" 2>&1 | tee "${ALIGN_LOG}"; then
    ALIGNMENT_STAGE_STATUS=0
  else
    ALIGNMENT_STAGE_STATUS=$?
  fi
fi
log "Stage alignment_pipeline (external) elapsed: $(($(date +%s) - ALIGNMENT_STAGE_TS))s" | tee -a "${ALIGN_LOG}"
if (( ALIGNMENT_STAGE_STATUS != 0 )); then exit "${ALIGNMENT_STAGE_STATUS}"; fi

ALIGNED_INDEX_TS=$(date +%s)
if samtools index "${ALIGNED_BAM}" 2>&1 | tee -a "${ALIGN_LOG}"; then
  ALIGNED_INDEX_STATUS=0
else
  ALIGNED_INDEX_STATUS=$?
fi
log "Stage aligned_bam_index (external) elapsed: $(($(date +%s) - ALIGNED_INDEX_TS))s" | tee -a "${ALIGN_LOG}"
if (( ALIGNED_INDEX_STATUS != 0 )); then exit "${ALIGNED_INDEX_STATUS}"; fi
step_end "Step 2/6"

step_start
log "Step 3/6: Sockeye-style dual-end barcode/UMI tags on aligned BAM"
pushd "${ALIGN_DIR}" >/dev/null

run_stage prepare_read_tags "${DOWNSTREAM_DIR}/prepare_read_tags.py" \
  --input "${READ_ASSIGNED_CELL}" \
  --output "${SAMPLE_ID}.read_tags.tsv" 2>&1 | tee "${TAG_LOG}"

if [[ "${FEATURE_BACKEND}" == "ont" || "${GENE_TAGGING_BACKEND}" != "rust" ]]; then
  run_stage add_cb_ur_tags "${DOWNSTREAM_DIR}/add_cb_ur_tags.py" \
    --bam "${SAMPLE_ID}.aligned.sorted.bam" \
    --tags "${SAMPLE_ID}.read_tags.tsv" \
    --output "${SAMPLE_ID}.filtered.cb_ur.sorted.bam" \
    --threads "${CLUSTER_THREADS}" 2>&1 | tee -a "${TAG_LOG}"

  if [[ "${FEATURE_BACKEND}" == "legacy" ]]; then
  BAMTOBED_TS=$(date +%s)
  if bedtools bamtobed -i "${SAMPLE_ID}.filtered.cb_ur.sorted.bam" > "${SAMPLE_ID}.filtered.cb_ur.bed"; then
    BAMTOBED_STATUS=0
  else
    BAMTOBED_STATUS=$?
  fi
  log "Stage bamtobed (external) elapsed: $(($(date +%s) - BAMTOBED_TS))s" | tee -a "${TAG_LOG}"
  if (( BAMTOBED_STATUS != 0 )); then exit "${BAMTOBED_STATUS}"; fi
  fi
else
  log "Fused Rust gene tagging selected; intermediate CB/UR BAM and BED are skipped" | tee -a "${TAG_LOG}"
fi

popd >/dev/null
step_end "Step 3/6"

step_start
log "Step 4/6: ${FEATURE_BACKEND} feature assignment, directional UMI clustering, and matrix generation"
pushd "${MATRIX_DIR}" >/dev/null

if [[ "${FEATURE_BACKEND}" == "ont" ]]; then
  ONT_ARGS=()
  if [[ -n "${ONT_TEMP_DIR}" ]]; then
    ONT_ARGS+=(--temp-dir "${ONT_TEMP_DIR}")
  fi
  if [[ "${SAVE_INTERMEDIATE}" -eq 1 || "${LIGHT_OUTPUT}" -eq 0 ]]; then
    ONT_ARGS+=(--save-intermediate)
  fi
  python3 "$(python_asset "${DOWNSTREAM_DIR}/ont_features.py")" \
    --bam "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam" \
    --tagged-bam "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam" \
    --genome "${GENE_FASTA}" --gtf "${ISOFORM_GTF}" \
    --output-dir "${OUT_DIR}/transcriptome" --sample "${SAMPLE_ID}" \
    --output-bam "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" \
    --gene-assigns "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" \
    --transcript-assigns "${SAMPLE_ID}.read_transcript_assigns.tsv" \
    --mapq "${GENE_ASSIGN_MAPQ}" --threads "${CLUSTER_THREADS}" --jobs "${ONT_JOBS}" \
    --stringtie-opts="${STRINGTIE_OPTS}" "${ONT_ARGS[@]}" 2>&1 | tee "${GENE_LOG}"
else
if [[ "${GENE_TAGGING_BACKEND}" != "rust" ]]; then
  run_stage assign_genes "${DOWNSTREAM_DIR}/assign_genes.py" \
    --output "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" \
    --mapq "${GENE_ASSIGN_MAPQ}" \
    --chunk_size "${GENE_ASSIGN_CHUNK_SIZE}" \
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.bed" \
    "${GENE_GTF}" 2>&1 | tee "${GENE_LOG}"

  run_stage add_gene_tags "${DOWNSTREAM_DIR}/add_gene_tags.py" \
    --output "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" \
    --threads "${CLUSTER_THREADS}" \
    "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam" \
    "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" 2>&1 | tee -a "${GENE_LOG}"

  GENE_INDEX_TS=$(date +%s)
  if samtools index -@ "${CLUSTER_THREADS}" "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam"; then
    GENE_INDEX_STATUS=0
  else
    GENE_INDEX_STATUS=$?
  fi
  log "Stage gene_bam_index (external) elapsed: $(($(date +%s) - GENE_INDEX_TS))s" | tee -a "${GENE_LOG}"
  if (( GENE_INDEX_STATUS != 0 )); then exit "${GENE_INDEX_STATUS}"; fi
fi

if [[ "${GENE_TAGGING_BACKEND}" == "rust" || "${GENE_TAGGING_BACKEND}" == "validate" ]]; then
  FUSED_ASSIGNMENTS="${SAMPLE_ID}.filtered.read_gene_assigns.tsv"
  FUSED_BAM="${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam"
  if [[ "${GENE_TAGGING_BACKEND}" == "validate" ]]; then
    FUSED_ASSIGNMENTS="${SAMPLE_ID}.rust.read_gene_assigns.tsv"
    FUSED_BAM="${SAMPLE_ID}.rust.cb_ur.gn.sorted.bam"
  fi
  run_stage tag_and_assign_genes "${DOWNSTREAM_DIR}/tag_and_assign_genes.py" \
    --bam "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam" \
    --tags "${ALIGN_DIR}/${SAMPLE_ID}.read_tags.tsv" \
    --gtf "${GENE_GTF}" \
    --gene-assigns "${FUSED_ASSIGNMENTS}" \
    --output "${FUSED_BAM}" \
    --mapq "${GENE_ASSIGN_MAPQ}" \
    --threads "${CLUSTER_THREADS}" 2>&1 | tee -a "${GENE_LOG}"
fi

if [[ "${GENE_TAGGING_BACKEND}" == "validate" ]]; then
  VALIDATE_DIR="$(mktemp -d "${MATRIX_DIR}/gene-tag-validation.XXXXXX")"
  cmp "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" "${SAMPLE_ID}.rust.read_gene_assigns.tsv"
  samtools view -h "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" > "${VALIDATE_DIR}/legacy.sam"
  samtools view -h "${SAMPLE_ID}.rust.cb_ur.gn.sorted.bam" > "${VALIDATE_DIR}/rust.sam"
  cmp "${VALIDATE_DIR}/legacy.sam" "${VALIDATE_DIR}/rust.sam"
  rm -f "${VALIDATE_DIR}/legacy.sam" "${VALIDATE_DIR}/rust.sam" \
    "${SAMPLE_ID}.rust.read_gene_assigns.tsv" "${SAMPLE_ID}.rust.cb_ur.gn.sorted.bam" \
    "${SAMPLE_ID}.rust.cb_ur.gn.sorted.bam.bai"
  rmdir "${VALIDATE_DIR}"
  log "Gene-tag validation passed: assignments and SAM records are identical" | tee -a "${GENE_LOG}"
fi
fi

run_stage cluster_umis_allbam "${DOWNSTREAM_DIR}/cluster_umis_allbam.py" \
  --output "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
  --ref_interval "${REF_INTERVAL}" \
  --cell_gene_max_reads "${CELL_GENE_MAX_READS}" \
  --threads "${CLUSTER_THREADS}" \
  "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" 2>&1 | tee -a "${GENE_LOG}"

run_stage cell_umi_gene_table "${DOWNSTREAM_DIR}/cell_umi_gene_table.py" \
  --output "${SAMPLE_ID}.cell_umi_gene.tsv" \
  --gene-expression-output "${SAMPLE_ID}.gene_expression.tsv" \
  "${SAMPLE_ID}.filtered.tagged.sorted.bam" 2>&1 | tee -a "${GENE_LOG}"

RNA_CLUSTER_TS=$(date +%s)
if python3 "$(python_asset "${DOWNSTREAM_DIR}/rna_cluster_analysis.py")" \
  --input "${SAMPLE_ID}.gene_expression.tsv" \
  --output "${SAMPLE_ID}.rna_cluster.tsv" 2>&1 | tee -a "${GENE_LOG}"; then
  RNA_CLUSTER_STATUS=0
else
  RNA_CLUSTER_STATUS=$?
fi
log "Stage rna_cluster_analysis (python) elapsed: $(($(date +%s) - RNA_CLUSTER_TS))s" | tee -a "${GENE_LOG}"
if (( RNA_CLUSTER_STATUS != 0 )); then exit "${RNA_CLUSTER_STATUS}"; fi

popd >/dev/null
step_end "Step 4/6"

step_start
if [[ "${SKIP_ISOFORM}" -eq 1 ]]; then
  log "Step 5/6: transcript assignment and isoform matrix generation skipped (--skip-isoform)"
  : > "${ISOFORM_LOG}"
else
  log "Step 5/6: transcript assignment and isoform matrix generation"
  pushd "${MATRIX_DIR}" >/dev/null

  if [[ "${FEATURE_BACKEND}" == "legacy" ]]; then
    run_stage assign_transcripts "${DOWNSTREAM_DIR}/assign_transcripts.py" \
    --output "${SAMPLE_ID}.read_transcript_assigns.tsv" \
    --mapq "${TRANSCRIPT_ASSIGN_MAPQ}" \
    --chunk_size "${TRANSCRIPT_ASSIGN_CHUNK_SIZE}" \
    "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
    "${ISOFORM_GTF}" 2>&1 | tee "${ISOFORM_LOG}"
  else
    log "Reusing ONT joint transcript assignments from Step 4/6" | tee "${ISOFORM_LOG}"
  fi

  run_stage isoform_expression "${DOWNSTREAM_DIR}/isoform_expression.py" \
    --output "${SAMPLE_ID}.isoform_expression.tsv" \
    "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
    "${SAMPLE_ID}.read_transcript_assigns.tsv" 2>&1 | tee -a "${ISOFORM_LOG}"

  popd >/dev/null
fi
step_end "Step 5/6"

if [[ "${FEATURE_BACKEND}" == "ont" ]]; then
  ONT_MATRIX_ARGS=()
  if [[ "${SKIP_ISOFORM}" -eq 0 ]]; then
    ONT_MATRIX_ARGS+=(--transcript-tsv "${MATRIX_DIR}/${SAMPLE_ID}.isoform_expression.tsv")
  fi
  python3 "$(python_asset "${DOWNSTREAM_DIR}/ont_matrices.py")" \
    --gene-tsv "${MATRIX_DIR}/${SAMPLE_ID}.gene_expression.tsv" \
    --output-dir "${MATRIX_DIR}" --sample "${SAMPLE_ID}" \
    "${ONT_MATRIX_ARGS[@]}" 2>&1 | tee -a "${ISOFORM_LOG}"
fi

step_start
log "Step 6/6: RNA QC and saturation analysis"
pushd "${QC_DIR}" >/dev/null

PARAMETERS_TSV="${SAMPLE_ID}.parameters.tsv"
write_parameters_tsv "${PARAMETERS_TSV}"

ln -sf "${MATRIX_DIR}/${SAMPLE_ID}.cell_umi_gene.tsv" cell_umi_gene.tsv
ln -sf "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam" filtered.sorted.bam
ln -sf "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam.bai" filtered.sorted.bam.bai

PRECOMPUTED_READ_QC_JSON="${UPSTREAM_DIR}/read_qc_summary.json"
PRECOMPUTED_FASTQ_COUNT="${UPSTREAM_DIR}/full_length_fastq_count.txt"
if [[ -s "${PRECOMPUTED_READ_QC_JSON}" && -s "${PRECOMPUTED_FASTQ_COUNT}" ]]; then
  cp "${PRECOMPUTED_READ_QC_JSON}" "${SAMPLE_ID}.read_qc_summary.json"
  cp "${PRECOMPUTED_FASTQ_COUNT}" "${SAMPLE_ID}.full_length_fastq_count.txt"
  log "Reused Read QC accumulated during the barcode FASTQ scan" | tee -a "${QC_LOG}"
else
  log "Precomputed Read QC unavailable; scanning the full-length FASTQ" | tee -a "${QC_LOG}"
  run_stage read_qc_summary "${DOWNSTREAM_DIR}/read_qc_summary.py" \
    --fastq "${FULL_LENGTH_FASTQ}" \
    --output-json "${SAMPLE_ID}.read_qc_summary.json" \
    --output-fastq-count "${SAMPLE_ID}.full_length_fastq_count.txt" 2>&1 | tee -a "${QC_LOG}"
fi

RAW_FASTQ_FOR_QC="${FULL_LENGTH_FASTQ}"
QC_THREADS=$(( THREADS < 8 ? THREADS : 8 ))

QC_ARGS=(
  --cell-umi-gene-tsv cell_umi_gene.tsv
  --bam filtered.sorted.bam
  --raw-fastq "${RAW_FASTQ_FOR_QC}"
  --full-length-fastq "${FULL_LENGTH_FASTQ}"
  --fastq-count-file "${SAMPLE_ID}.full_length_fastq_count.txt"
  --threads "${QC_THREADS}"
  --mixed-species
  --read-tags "${ALIGN_DIR}/${SAMPLE_ID}.read_tags.tsv"
  --barcode-validity-tsv "${BARCODE_VALIDITY_SUMMARY}"
)
if [[ "${SKIP_ISOFORM}" -eq 0 ]]; then
  QC_ARGS+=(--transcript-assigns "${MATRIX_DIR}/${SAMPLE_ID}.read_transcript_assigns.tsv")
fi
if [[ "${SKIP_GLYCINE}" -eq 0 ]]; then
  QC_ARGS+=(--glycine-log "${LOG_DIR}/00_glycine.log" --glycine-stats "${GLYCINE_STATS}")
fi

run_stage rna_qc_metrics "${DOWNSTREAM_DIR}/rna_qc_metrics_mixed.py" \
  "${QC_ARGS[@]}" 2>&1 | tee "${QC_LOG}"
if [[ ! -f rna_violin_plot.png ]]; then
  python3 "$(python_asset "${DOWNSTREAM_DIR}/rna_violin_plot.py")" \
    --input per_cell_qc.tsv \
    --output rna_violin_plot.png 2>&1 | tee -a "${QC_LOG}"
fi
mv -f rna_qc_metrics.tsv "${SAMPLE_ID}.rna_qc_metrics.tsv"
mv -f single_cell_report_metrics.tsv "${SAMPLE_ID}.single_cell_report_metrics.tsv"
mv -f per_cell_qc.tsv "${SAMPLE_ID}.per_cell_qc.tsv"
mv -f rna_violin_plot.png "${SAMPLE_ID}.rna_violin_plot.png"

run_stage saturation "${DOWNSTREAM_DIR}/Saturation.py" \
  --input cell_umi_gene.tsv \
  --output-tsv "${SAMPLE_ID}.saturation.tsv" \
  --output-png "${SAMPLE_ID}.saturation_curves.png" 2>&1 | tee -a "${QC_LOG}"
if [[ ! -f "${SAMPLE_ID}.saturation_curves.png" ]]; then
  python3 "$(python_asset "${DOWNSTREAM_DIR}/Saturation.py")" \
    --plot-existing-tsv "${SAMPLE_ID}.saturation.tsv" \
    --output-png "${SAMPLE_ID}.saturation_curves.png" 2>&1 | tee -a "${QC_LOG}"
fi

METRICS_SUMMARY_ARGS=(
  --sample-id "${SAMPLE_ID}"
  --species "$(basename "${REF_DIR:-$(dirname "${GENE_FASTA}")}")"
  --read-qc-json "${SAMPLE_ID}.read_qc_summary.json"
  --rna-qc-metrics-tsv "${SAMPLE_ID}.rna_qc_metrics.tsv"
  --saturation-tsv "${SAMPLE_ID}.saturation.tsv"
  --output-xlsx "metrics_summary.xlsx"
)
if [[ "${SKIP_GLYCINE}" -eq 1 ]]; then
  METRICS_SUMMARY_ARGS+=(--skip-glycine)
else
  METRICS_SUMMARY_ARGS+=(--glycine-stats "${GLYCINE_STATS}")
fi
run_stage metrics_summary "${DOWNSTREAM_DIR}/metrics_summary.py" \
  "${METRICS_SUMMARY_ARGS[@]}" 2>&1 | tee -a "${QC_LOG}"

run_stage barnyard_qc "${DOWNSTREAM_DIR}/barnyard_qc.py" \
  --input cell_umi_gene.tsv \
  --out-dir barnyard_qc \
  --ambient-umi-threshold 1000 \
  --singlet-threshold "${SINGLET_THRESHOLD}" 2>&1 | tee -a "${QC_LOG}"

BARNYARD_SUMMARY_TSV="barnyard_qc/barnyard_summary.tsv"
BARNYARD_PER_CELL_TSV="barnyard_qc/barnyard_per_cell.tsv"

BUILD_REPORT_ARGS=(
  --sample-id "${SAMPLE_ID}"
  --output-html "${SAMPLE_ID}.single_cell_report.html"
  --report-metrics-tsv "${SAMPLE_ID}.single_cell_report_metrics.tsv"
  --rna-qc-metrics-tsv "${SAMPLE_ID}.rna_qc_metrics.tsv"
  --saturation-tsv "${SAMPLE_ID}.saturation.tsv"
  --read-qc-json "${SAMPLE_ID}.read_qc_summary.json"
  --parameters-tsv "${PARAMETERS_TSV}"
  --per-cell-qc-tsv "${SAMPLE_ID}.per_cell_qc.tsv"
  --rna-cluster-tsv "${MATRIX_DIR}/${SAMPLE_ID}.rna_cluster.tsv"
  --barcode-counts-3p-tsv "${UPSTREAM_DIR}/barcode_counts_3p.tsv"
  --barcode-counts-5p-tsv "${UPSTREAM_DIR}/barcode_counts_5p.tsv"
  --whitelist-3p "${UPSTREAM_DIR}/whitelist_3p.csv"
  --whitelist-5p "${UPSTREAM_DIR}/whitelist_5p.csv"
  --read-assigned-cell "${READ_ASSIGNED_CELL}"
  --barcode-to-cell "${UPSTREAM_DIR}/barcode_to_cell.csv"
  --knee-plot-3p "${UPSTREAM_DIR}/knee_plot_3p.png"
  --knee-plot-5p "${UPSTREAM_DIR}/knee_plot_5p.png"
  --saturation-png "${SAMPLE_ID}.saturation_curves.png"
  --rna-violin-png "${SAMPLE_ID}.rna_violin_plot.png"
)
if [[ -f "${BARNYARD_SUMMARY_TSV}" && -f "${BARNYARD_PER_CELL_TSV}" ]]; then
  BUILD_REPORT_ARGS+=(
    --barnyard-summary-tsv "${BARNYARD_SUMMARY_TSV}"
    --barnyard-per-cell-tsv "${BARNYARD_PER_CELL_TSV}"
  )
else
  log "WARNING: Barnyard report inputs are incomplete; omitting Barnyard QC from the HTML report."
fi
if [[ "${SKIP_GLYCINE}" -eq 1 ]]; then
  BUILD_REPORT_ARGS+=(--skip-glycine)
else
  BUILD_REPORT_ARGS+=(--glycine-stats "${GLYCINE_STATS}")
fi
run_stage build_report "${DOWNSTREAM_DIR}/build_report.py" \
  "${BUILD_REPORT_ARGS[@]}" 2>&1 | tee -a "${QC_LOG}"

popd >/dev/null
step_end "Step 6/6"

if [[ "${LIGHT_OUTPUT}" -eq 1 && "${SAVE_INTERMEDIATE}" -eq 0 ]]; then
  log "Light output: removing reproducible BAM and read-level intermediates"
  cleanup_large_intermediates
fi
if [[ "${REMOVE_FINAL_BAM}" -eq 1 ]]; then
  log "--remove-final-bam: removing all alignment and matrix BAM outputs"
  cleanup_all_bam_outputs
fi

log "Pipeline completed successfully."
SCRIPT_END_TS=$(date +%s)
log "Total elapsed: $((SCRIPT_END_TS - SCRIPT_START_TS))s"
log "Key outputs:"
echo "  full-length fastq        : ${FULL_LENGTH_FASTQ}"
echo "  barcode-to-cell mapping  : ${UPSTREAM_DIR}/barcode_to_cell.csv"
if [[ "${SKIP_CELL_FASTQ}" -eq 0 ]]; then
  echo "  cell reads fastq         : ${CELL_READS_FASTQ}"
fi
if [[ "${SKIP_MATCHED_FASTQ}" -eq 0 ]]; then
  echo "  matched reads fastq      : ${MATCHED_READS_FASTQ}"
fi
if [[ "${REMOVE_FINAL_BAM}" -eq 0 ]]; then
  echo "  tagged bam              : ${MATRIX_DIR}/${SAMPLE_ID}.filtered.tagged.sorted.bam"
else
  echo "  tagged bam              : removed (--remove-final-bam)"
fi
echo "  gene expression         : ${MATRIX_DIR}/${SAMPLE_ID}.gene_expression.tsv"
echo "  RNA cluster table       : ${MATRIX_DIR}/${SAMPLE_ID}.rna_cluster.tsv"
if [[ "${SKIP_ISOFORM}" -eq 0 ]]; then
  echo "  isoform expression      : ${MATRIX_DIR}/${SAMPLE_ID}.isoform_expression.tsv"
else
  echo "  isoform expression      : skipped (--skip-isoform)"
fi
echo "  RNA QC                  : ${QC_DIR}/${SAMPLE_ID}.rna_qc_metrics.tsv"
echo "  metrics workbook        : ${QC_DIR}/metrics_summary.xlsx"
echo "  report metrics         : ${QC_DIR}/${SAMPLE_ID}.single_cell_report_metrics.tsv"
echo "  html report            : ${QC_DIR}/${SAMPLE_ID}.single_cell_report.html"
echo "  saturation table        : ${QC_DIR}/${SAMPLE_ID}.saturation.tsv"
