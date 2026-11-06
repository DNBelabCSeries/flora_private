use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use tempfile::TempDir;

#[derive(Clone, Copy, Debug)]
pub enum WorkflowKind {
    Single,
    Mixed,
}

const SINGLE_RUNNER: &[u8] = include_bytes!("../run_all.sh");
const MIXED_RUNNER: &[u8] = include_bytes!("../run_all_mixed_species.sh");

macro_rules! pyc {
    ($name:literal) => {
        include_bytes!(concat!(env!("OUT_DIR"), "/flora-python/", $name))
    };
}

const PYTHON_ASSETS: &[(&str, &[u8])] = &[
    ("Saturation.pyc", pyc!("Saturation.pyc")),
    ("add_cb_ur_tags.pyc", pyc!("add_cb_ur_tags.pyc")),
    ("add_gene_tags.pyc", pyc!("add_gene_tags.pyc")),
    ("assign_genes.pyc", pyc!("assign_genes.pyc")),
    ("assign_transcripts.pyc", pyc!("assign_transcripts.pyc")),
    ("barnyard_qc.pyc", pyc!("barnyard_qc.pyc")),
    (
        "build_mixed_species_reference.pyc",
        pyc!("build_mixed_species_reference.pyc"),
    ),
    ("build_report.pyc", pyc!("build_report.pyc")),
    ("cell_umi_gene_table.pyc", pyc!("cell_umi_gene_table.pyc")),
    ("cluster_umis_allbam.pyc", pyc!("cluster_umis_allbam.pyc")),
    ("gene_expression.pyc", pyc!("gene_expression.pyc")),
    (
        "generate_26bp_whitelists.pyc",
        pyc!("generate_26bp_whitelists.pyc"),
    ),
    ("generate_knee_plots.pyc", pyc!("generate_knee_plots.pyc")),
    ("isoform_expression.pyc", pyc!("isoform_expression.pyc")),
    ("metrics_summary.pyc", pyc!("metrics_summary.pyc")),
    ("ont_features.pyc", pyc!("ont_features.pyc")),
    ("ont_matrices.pyc", pyc!("ont_matrices.pyc")),
    ("prepare_read_tags.pyc", pyc!("prepare_read_tags.pyc")),
    ("read_qc_summary.pyc", pyc!("read_qc_summary.pyc")),
    ("rna_cluster_analysis.pyc", pyc!("rna_cluster_analysis.pyc")),
    ("rna_qc_metrics.pyc", pyc!("rna_qc_metrics.pyc")),
    ("rna_qc_metrics_mixed.pyc", pyc!("rna_qc_metrics_mixed.pyc")),
    ("rna_violin_plot.pyc", pyc!("rna_violin_plot.pyc")),
];

const STATIC_ASSETS: &[(&str, &[u8])] = &[
    (
        "report_template.html",
        include_bytes!("../scripts/report_template.html"),
    ),
    (
        "plotly-2.26.0.min.js",
        include_bytes!("../scripts/plotly-2.26.0.min.js"),
    ),
];

const INTERNAL_STAGE_NAMES: &[&str] = &[
    "add_cb_ur_tags",
    "add_gene_tags",
    "assign_genes",
    "assign_transcripts",
    "cell_umi_gene_table",
    "cluster_umis_allbam",
    "gene_expression",
    "generate_26bp_whitelists",
    "isoform_expression",
    "prepare_read_tags",
    "read_qc_summary",
    "rna_qc_metrics",
    "saturation",
    "tag_and_assign_genes",
];

pub fn run(kind: WorkflowKind, mut args: Vec<OsString>) -> Result<()> {
    let explicit_python = take_option(&mut args, "--python")?;
    let python = resolve_python(explicit_python)?;
    normalize_path_options(&mut args)?;
    let runtime = prepare_runtime(kind)?;
    let runner = runtime.path().join(match kind {
        WorkflowKind::Single => "run_all.sh",
        WorkflowKind::Mixed => "run_all_mixed_species.sh",
    });

    let old_path = std::env::var_os("PATH").unwrap_or_default();
    let joined_path = std::env::join_paths(
        std::iter::once(runtime.path().join("bin")).chain(std::env::split_paths(&old_path)),
    )?;
    let status = Command::new("bash")
        .arg(&runner)
        .args(args)
        .env("PATH", joined_path)
        .env("FLORA_PYTHON", &python)
        .env("FLORA_INTERNAL_WORKFLOW", "1")
        .env("MPLCONFIGDIR", runtime.path().join("matplotlib"))
        .status()
        .with_context(|| format!("failed to start embedded workflow {}", runner.display()))?;

    if !status.success() {
        match status.code() {
            Some(code) => bail!("Flora workflow failed with exit code {code}"),
            None => bail!("Flora workflow was terminated by a signal"),
        }
    }
    Ok(())
}

fn prepare_runtime(kind: WorkflowKind) -> Result<TempDir> {
    let runtime = tempfile::Builder::new()
        .prefix("flora-runtime-")
        .tempdir()
        .context("failed to create Flora private runtime directory")?;
    let root = runtime.path();
    let scripts = root.join("scripts");
    let release = root.join("target/release");
    let bin = root.join("bin");
    fs::create_dir_all(&scripts)?;
    fs::create_dir_all(&release)?;
    fs::create_dir_all(&bin)?;
    fs::create_dir_all(root.join("matplotlib"))?;

    write_runner(&root.join("run_all.sh"), SINGLE_RUNNER, false)?;
    write_runner(&root.join("run_all_mixed_species.sh"), MIXED_RUNNER, true)?;
    write_asset(&root.join("main.pyc"), pyc!("main.pyc"))?;
    for (name, content) in PYTHON_ASSETS {
        write_asset(&scripts.join(name), content)?;
    }
    for (name, content) in STATIC_ASSETS {
        write_asset(&scripts.join(name), content)?;
    }

    let executable = std::env::current_exe().context("failed to locate the Flora executable")?;
    make_link(&executable, &release.join("flora"))?;
    for stage in INTERNAL_STAGE_NAMES {
        make_link(&executable, &release.join(stage))?;
    }

    let wrapper = bin.join("python3");
    write_asset(
        &wrapper,
        b"#!/usr/bin/env bash\nexec \"$FLORA_PYTHON\" \"$@\"\n",
    )?;
    make_executable(&wrapper)?;

    let selected_runner = match kind {
        WorkflowKind::Single => root.join("run_all.sh"),
        WorkflowKind::Mixed => root.join("run_all_mixed_species.sh"),
    };
    if !selected_runner.is_file() {
        bail!("embedded workflow runner is missing")
    }
    Ok(runtime)
}

fn write_asset(path: &Path, content: &[u8]) -> Result<()> {
    fs::write(path, content).with_context(|| format!("failed to extract {}", path.display()))
}

fn write_runner(path: &Path, content: &[u8], mixed: bool) -> Result<()> {
    let mut text = String::from_utf8(content.to_vec()).context("embedded runner is not UTF-8")?;
    text = if mixed {
        text.replace("bash run_all_mixed_species.sh", "flora mixed")
            .replace("bash run_all.sh", "flora mixed")
    } else {
        text.replace("bash run_all.sh", "flora")
    };
    write_asset(path, text.as_bytes())
}

#[cfg(unix)]
fn make_link(target: &Path, link: &Path) -> Result<()> {
    std::os::unix::fs::symlink(target, link)
        .with_context(|| format!("failed to link internal stage {}", link.display()))
}

#[cfg(not(unix))]
fn make_link(_target: &Path, _link: &Path) -> Result<()> {
    bail!("Flora full workflows currently require Linux or macOS")
}

#[cfg(unix)]
fn make_executable(path: &Path) -> Result<()> {
    use std::os::unix::fs::PermissionsExt;
    let mut permissions = fs::metadata(path)?.permissions();
    permissions.set_mode(0o700);
    fs::set_permissions(path, permissions)?;
    Ok(())
}

#[cfg(not(unix))]
fn make_executable(_path: &Path) -> Result<()> {
    Ok(())
}

fn take_option(args: &mut Vec<OsString>, name: &str) -> Result<Option<OsString>> {
    let mut selected = None;
    let mut index = 0;
    while index < args.len() {
        let value = args[index].to_string_lossy();
        if value == name {
            if index + 1 >= args.len() {
                bail!("{name} requires a value")
            }
            selected = Some(args.remove(index + 1));
            args.remove(index);
            continue;
        }
        if let Some(value) = value.strip_prefix(&format!("{name}=")) {
            selected = Some(OsString::from(value));
            args.remove(index);
            continue;
        }
        index += 1;
    }
    Ok(selected)
}

fn resolve_python(explicit: Option<OsString>) -> Result<PathBuf> {
    let requested = explicit
        .or_else(|| std::env::var_os("FLORA_PYTHON"))
        .unwrap_or_else(|| OsString::from("python3"));
    let requested_path = PathBuf::from(&requested);
    if requested_path.components().count() > 1 {
        if requested_path.is_file() {
            return Ok(requested_path);
        }
        bail!("Python interpreter not found: {}", requested_path.display())
    }
    let path = std::env::var_os("PATH").unwrap_or_default();
    for directory in std::env::split_paths(&path) {
        let candidate = directory.join(&requested_path);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }
    bail!(
        "Python interpreter not found in PATH: {}",
        requested_path.display()
    )
}

fn normalize_path_options(args: &mut [OsString]) -> Result<()> {
    const PATH_OPTIONS: &[&str] = &[
        "--fastq-dir",
        "--full-length-fastq",
        "--barcode-list-10bp",
        "--barcode_list_10bp",
        "--ref-dir",
        "--ref_dir",
        "--gene-fasta",
        "--genome-fa",
        "--junction-bed",
        "--chrom-sizes",
        "--gene-gtf",
        "--isoform-gtf",
        "--out-dir",
        "--outdir",
        "--glycine-outdir",
        "--ont-temp-dir",
    ];
    let cwd = std::env::current_dir().context("failed to resolve current directory")?;
    let mut index = 0;
    while index < args.len() {
        let raw = args[index].to_string_lossy();
        if raw == "--fastq" {
            index += 1;
            let first_value = index;
            while index < args.len() && !args[index].to_string_lossy().starts_with('-') {
                let path = PathBuf::from(&args[index]);
                if path.is_relative() {
                    args[index] = cwd.join(path).into_os_string();
                }
                index += 1;
            }
            if index == first_value {
                bail!("--fastq requires at least one value")
            }
            continue;
        }
        if PATH_OPTIONS.iter().any(|option| raw == *option) {
            if index + 1 >= args.len() {
                bail!("{} requires a value", raw)
            }
            let path = PathBuf::from(&args[index + 1]);
            if path.is_relative() {
                args[index + 1] = cwd.join(path).into_os_string();
            }
            index += 2;
            continue;
        }
        if let Some((option, value)) = raw.split_once('=') {
            if option == "--fastq" || PATH_OPTIONS.contains(&option) {
                let path = PathBuf::from(value);
                if path.is_relative() {
                    args[index] = OsString::from(format!("{option}={}", cwd.join(path).display()));
                }
            }
        }
        index += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{normalize_path_options, MIXED_RUNNER, SINGLE_RUNNER};
    use std::ffi::OsString;
    use std::path::PathBuf;

    #[test]
    fn normalizes_every_multi_fastq_value_and_fastq_directory() {
        let cwd = std::env::current_dir().unwrap();
        let mut files = vec![
            OsString::from("--fastq"),
            OsString::from("sample_1.fastq.gz"),
            OsString::from("sample_2.fastq.gz"),
            OsString::from("--out-dir"),
            OsString::from("output"),
        ];
        normalize_path_options(&mut files).unwrap();
        assert_eq!(PathBuf::from(&files[1]), cwd.join("sample_1.fastq.gz"));
        assert_eq!(PathBuf::from(&files[2]), cwd.join("sample_2.fastq.gz"));

        let mut directory = vec![OsString::from("--fastq-dir"), OsString::from("chunks")];
        normalize_path_options(&mut directory).unwrap();
        assert_eq!(PathBuf::from(&directory[1]), cwd.join("chunks"));
    }

    #[test]
    fn embedded_workflows_default_annotation_mapq_to_30() {
        for runner in [SINGLE_RUNNER, MIXED_RUNNER] {
            let runner = std::str::from_utf8(runner).unwrap();
            assert!(runner.contains("GENE_ASSIGN_MAPQ=30"));
            assert!(runner.contains("TRANSCRIPT_ASSIGN_MAPQ=30"));
        }
    }
}
