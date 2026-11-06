use std::ffi::OsString;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use clap::{ArgAction, Parser};
use flora::pipeline::{run_pipeline, PipelineConfig};
use flora::workflow_runtime::{self, WorkflowKind};

#[derive(Debug, Parser)]
#[command(version, about = "Barcode, UMI, and cell assignment for Flora")]
struct Cli {
    #[arg(required = true)]
    fastq_fns: Vec<PathBuf>,

    #[arg(long = "full-bc-whitelist-3p")]
    full_bc_whitelist_3p: PathBuf,

    #[arg(long = "full-bc-whitelist-5p")]
    full_bc_whitelist_5p: PathBuf,

    #[arg(long = "no-revcomp-whitelist", action = ArgAction::SetTrue)]
    no_revcomp_whitelist: bool,

    #[arg(long = "out_dir", default_value = "Flora")]
    out_dir: PathBuf,

    #[arg(long = "threads", default_value_t = 32)]
    threads: usize,

    #[arg(long = "batch_size", default_value_t = 100_000)]
    batch_size: usize,

    #[arg(long = "assign_batchsize", default_value_t = 10_000)]
    assign_batchsize: usize,

    #[arg(long = "minQ", default_value_t = 2)]
    min_q: i32,

    #[arg(long = "exp_cells", default_value_t = 5000)]
    exp_cells: usize,

    #[arg(long = "max_ed", default_value_t = 2)]
    max_ed: usize,

    #[arg(long = "PAIR_MIN")]
    pair_min: Option<usize>,

    #[arg(long = "auto_pair_min_floor", default_value_t = 10)]
    auto_pair_min_floor: usize,

    #[arg(long = "auto_pair_min_quantile", default_value_t = 0.1)]
    auto_pair_min_quantile: f64,

    #[arg(long = "dominance_min", default_value_t = 0.80)]
    dominance_min: f64,

    #[arg(long = "drop_umiA_ratio_gt", default_value_t = 0.5)]
    drop_umi_a_ratio_gt: f64,

    #[arg(long = "TOP1_ALPHA", default_value_t = 0.1)]
    top1_alpha: f64,

    #[arg(long = "TOP1_ALPHA_UMI", default_value_t = 0.3)]
    top1_alpha_umi: f64,

    #[arg(long = "min_reads_per_cell", default_value_t = 20)]
    min_reads_per_cell: usize,

    #[arg(long = "BC_fixed_3p", default_value = "GGTAGC")]
    bc_fixed_3p: String,

    #[arg(long = "umi_fixed_3p", default_value = "GATCT")]
    umi_fixed_3p: String,

    #[arg(long = "BC_fixed_5p", default_value = "GGAAGG")]
    bc_fixed_5p: String,

    #[arg(long = "umi_fixed_5p", default_value = "CAGCA")]
    umi_fixed_5p: String,

    #[arg(long = "barcode_extract_mode", default_value = "fixed_seq")]
    barcode_extract_mode: String,

    #[arg(long = "light-output", action = ArgAction::SetTrue, default_value_t = true)]
    light_output: bool,

    #[arg(long = "full-output", action = ArgAction::SetTrue)]
    full_output: bool,

    #[arg(long = "alignment-fastq", hide = true)]
    alignment_fastq: Option<PathBuf>,

    #[arg(long = "save-intermediate", action = ArgAction::SetTrue)]
    save_intermediate: bool,

    #[arg(long = "save_merge_debug", action = ArgAction::SetTrue)]
    save_merge_debug: bool,

    #[arg(long = "require_pass_both_ends", action = ArgAction::SetTrue)]
    require_pass_both_ends: bool,

    #[arg(long = "include_other_components", action = ArgAction::SetTrue, default_value_t = false)]
    include_other_components: bool,

    #[arg(long = "exclude_other_components", action = ArgAction::SetTrue)]
    exclude_other_components: bool,

    #[arg(long = "max_other_component_barcodes", default_value_t = 8)]
    max_other_component_barcodes: usize,

    #[arg(long = "absorb_unassigned_paired", action = ArgAction::SetTrue, default_value_t = true)]
    absorb_unassigned_paired: bool,

    #[arg(long = "disable_absorb_unassigned_paired", action = ArgAction::SetTrue)]
    disable_absorb_unassigned_paired: bool,
}

fn main() -> Result<()> {
    let mut args: Vec<OsString> = std::env::args_os().collect();
    let executable_name = args
        .first()
        .and_then(|value| std::path::Path::new(value).file_name())
        .and_then(|value| value.to_str())
        .unwrap_or("flora");
    if let Some(result) = dispatch_internal_stage(executable_name) {
        return result;
    }

    match args.get(1).and_then(|value| value.to_str()) {
        Some("glycine") => {
            args.remove(1);
            args[0] = OsString::from("flora glycine");
            return flora::glycine::batch::run_from(args);
        }
        Some("analyze") => {
            args.remove(1);
            args[0] = OsString::from("flora analyze");
            return run_analyze_from(args);
        }
        Some("mixed") => {
            args.drain(0..2);
            return workflow_runtime::run(WorkflowKind::Mixed, args);
        }
        Some("run") | Some("run-mixed") => {
            anyhow::bail!(
                "unsupported command; use 'flora' for single species or 'flora mixed' for mixed species"
            );
        }
        Some("help") if args.len() == 2 => {
            print_top_level_help();
            return workflow_runtime::run(WorkflowKind::Single, vec![OsString::from("--help")]);
        }
        Some("--help") | Some("-h") if args.len() == 2 => {
            print_top_level_help();
            return workflow_runtime::run(WorkflowKind::Single, vec![OsString::from("--help")]);
        }
        Some("--version") | Some("-V") if args.len() == 2 => {
            println!("flora {}", env!("CARGO_PKG_VERSION"));
            return Ok(());
        }
        _ if std::env::var_os("FLORA_INTERNAL_WORKFLOW").is_some() => run_analyze_from(args),
        _ => {
            args.remove(0);
            workflow_runtime::run(WorkflowKind::Single, args)
        }
    }
}

fn print_top_level_help() {
    println!("Flora: end-to-end full-length single-cell transcriptomics\n");
    println!("Commands:");
    println!("  flora [OPTIONS]          Run the complete single-species workflow");
    println!("  flora mixed [OPTIONS]    Run the complete mixed-species workflow");
    println!("  flora glycine [OPTIONS]  Run Glycine full-length read identification only");
    println!("  flora analyze [OPTIONS]  Run barcode, UMI, and cell assignment only\n");
}

fn dispatch_internal_stage(name: &str) -> Option<Result<()>> {
    match name {
        "add_cb_ur_tags" => Some(flora::stage_add_cb_ur_tags::main()),
        "add_gene_tags" => Some(flora::stage_add_gene_tags::main()),
        "assign_genes" => Some(flora::stage_assign_genes::main()),
        "assign_transcripts" => Some(flora::stage_assign_transcripts::main()),
        "cell_umi_gene_table" => Some(flora::stage_cell_umi_gene_table::main()),
        "cluster_umis_allbam" => Some(flora::stage_cluster_umis_allbam::main()),
        "gene_expression" => Some(flora::stage_gene_expression::main()),
        "generate_26bp_whitelists" => Some(flora::stage_generate_26bp_whitelists::main()),
        "isoform_expression" => Some(flora::stage_isoform_expression::main()),
        "prepare_read_tags" => Some(flora::stage_prepare_read_tags::main()),
        "read_qc_summary" => Some(flora::stage_read_qc_summary::main()),
        "rna_qc_metrics" => Some(flora::stage_rna_qc_metrics::main()),
        "saturation" => Some(flora::stage_saturation::main()),
        "tag_and_assign_genes" => Some(flora::stage_tag_and_assign_genes::main()),
        _ => None,
    }
}

fn run_analyze_from<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let total_t0 = Instant::now();
    let cli = Cli::parse_from(args);
    if cli.barcode_extract_mode != "fixed_seq" {
        anyhow::bail!("Flora currently supports --barcode_extract_mode fixed_seq only");
    }
    if cli.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global();
    }

    let skip_fastq = cli.light_output && !cli.full_output;
    let include_other_components = cli.include_other_components && !cli.exclude_other_components;
    let absorb_unassigned_paired =
        cli.absorb_unassigned_paired && !cli.disable_absorb_unassigned_paired;
    let config = PipelineConfig {
        fastq_fns: cli.fastq_fns,
        full_bc_whitelist_3p: cli.full_bc_whitelist_3p,
        full_bc_whitelist_5p: cli.full_bc_whitelist_5p,
        revcomp_whitelist: !cli.no_revcomp_whitelist,
        out_dir: cli.out_dir,
        batch_size: cli.batch_size,
        assign_batchsize: cli.assign_batchsize,
        min_q: cli.min_q,
        exp_cells: cli.exp_cells,
        max_ed: cli.max_ed,
        pair_min: cli.pair_min,
        auto_pair_min_floor: cli.auto_pair_min_floor,
        auto_pair_min_quantile: cli.auto_pair_min_quantile,
        dominance_min: cli.dominance_min,
        drop_umi_a_ratio_gt: cli.drop_umi_a_ratio_gt,
        top1_alpha: cli.top1_alpha,
        top1_alpha_umi: cli.top1_alpha_umi,
        require_pass_both_ends: cli.require_pass_both_ends,
        include_other_components,
        max_other_component_barcodes: cli.max_other_component_barcodes,
        absorb_unassigned_paired,
        min_reads_per_cell: cli.min_reads_per_cell,
        bc_fixed_3p: cli.bc_fixed_3p,
        umi_fixed_3p: cli.umi_fixed_3p,
        bc_fixed_5p: cli.bc_fixed_5p,
        umi_fixed_5p: cli.umi_fixed_5p,
        skip_matched_fastq: skip_fastq,
        skip_unmatched_fastq: skip_fastq,
        skip_cell_fastq: skip_fastq,
        alignment_fastq: cli.alignment_fastq,
        save_intermediate: cli.save_intermediate,
        save_correction_maps: cli.save_merge_debug,
    };

    let summary = run_pipeline(&config)?;
    print_summary(&summary, &config.out_dir);
    println!(
        "\n[Flora] Total elapsed: {:.2}s",
        total_t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn print_summary(summary: &flora::pipeline::PipelineSummary, out_dir: &std::path::Path) {
    println!("\n=== Flora Summary ===");
    println!("Output directory: {}", out_dir.display());
    println!("FASTQ files: {}", summary.fastq_files);
    println!("Reads total: {}", summary.reads_total);
    println!("Reads demultiplexed: {}", summary.reads_demultiplexed);
    println!("\nBarcode validity in corrected reads:");
    println!(
        "Any barcode valid: {} ({:.2}%)",
        summary.barcode_validity_stats.valid_any_n,
        summary.barcode_validity_stats.valid_any_ratio * 100.0
    );
    println!(
        "Only 3' valid: {} ({:.2}%)",
        summary.barcode_validity_stats.only_3p_n,
        summary.barcode_validity_stats.only_3p_ratio * 100.0
    );
    println!(
        "Only 5' valid: {} ({:.2}%)",
        summary.barcode_validity_stats.only_5p_n,
        summary.barcode_validity_stats.only_5p_ratio * 100.0
    );
    println!(
        "Both ends valid: {} ({:.2}%)",
        summary.barcode_validity_stats.both_n,
        summary.barcode_validity_stats.both_ratio * 100.0
    );
    println!(
        "Neither valid: {} ({:.2}%)",
        summary.barcode_validity_stats.neither_n,
        summary.barcode_validity_stats.neither_ratio * 100.0
    );
    println!("Merged putative rows: {}", summary.merged_rows);
    println!("Filtered rows: {}", summary.filtered_rows);
    println!(
        "Unique BC3_20bp_rc: {}",
        summary.trimmed_barcode_uniques.unique_bc3_20bp_rc
    );
    println!(
        "Unique BC5_20bp: {}",
        summary.trimmed_barcode_uniques.unique_bc5_20bp
    );
    println!(
        "Unique union (BC3 ∪ BC5): {}",
        summary.trimmed_barcode_uniques.unique_union_bc3_bc5
    );
    println!("Clean read rows: {}", summary.clean_reads_rows);
    println!("Pair rows kept: {}", summary.pair_counts_rows);
    println!(
        "PAIR_MIN used: {} (mode={})",
        summary.pair_stats.pair_min, summary.pair_stats.pair_min_mode
    );
    println!(
        "TOP1 alpha (reads/UMI): {} / {}",
        summary.pair_stats.top1_alpha, summary.pair_stats.top1_alpha_umi
    );
    println!(
        "Core cells: {} core barcodes: {}",
        summary.core_cells, summary.core_barcodes
    );
    println!("Assigned read rows: {}", summary.assigned_rows);
    println!("\nassign_stats:");
    println!("{:?}", summary.assign_stats);
    println!("\nTop cell_read_stats:");
    for row in &summary.cell_read_stats_head {
        println!("{:?}", row);
    }
}
