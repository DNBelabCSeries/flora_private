use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const PYTHON_SOURCES: &[&str] = &[
    "main.py",
    "scripts/Saturation.py",
    "scripts/add_cb_ur_tags.py",
    "scripts/add_gene_tags.py",
    "scripts/assign_genes.py",
    "scripts/assign_transcripts.py",
    "scripts/barnyard_qc.py",
    "scripts/build_mixed_species_reference.py",
    "scripts/build_report.py",
    "scripts/cell_umi_gene_table.py",
    "scripts/cluster_umis_allbam.py",
    "scripts/gene_expression.py",
    "scripts/generate_26bp_whitelists.py",
    "scripts/generate_knee_plots.py",
    "scripts/isoform_expression.py",
    "scripts/metrics_summary.py",
    "scripts/ont_features.py",
    "scripts/ont_matrices.py",
    "scripts/prepare_read_tags.py",
    "scripts/read_qc_summary.py",
    "scripts/rna_cluster_analysis.py",
    "scripts/rna_qc_metrics.py",
    "scripts/rna_qc_metrics_mixed.py",
    "scripts/rna_violin_plot.py",
];

fn main() {
    println!("cargo:rerun-if-env-changed=PYTHON_BIN");
    let root = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("flora-python");
    fs::create_dir_all(&output).expect("failed to create embedded Python output directory");
    let python = env::var_os("PYTHON_BIN").unwrap_or_else(|| "python3".into());

    for relative in PYTHON_SOURCES {
        let source = root.join(relative);
        println!("cargo:rerun-if-changed={}", source.display());
        let file_name = Path::new(relative).file_name().unwrap().to_string_lossy();
        let destination = output.join(format!("{file_name}c"));
        let display_name = format!("flora-runtime/{relative}");
        let code = r#"import py_compile, sys
if sys.version_info[:2] != (3, 11):
    raise SystemExit(f"Flora embedded assets require Python 3.11, got {sys.version.split()[0]}")
py_compile.compile(
    sys.argv[1],
    cfile=sys.argv[2],
    dfile=sys.argv[3],
    doraise=True,
    invalidation_mode=py_compile.PycInvalidationMode.UNCHECKED_HASH,
)
"#;
        let status = Command::new(&python)
            .args(["-c", code])
            .arg(&source)
            .arg(&destination)
            .arg(&display_name)
            .status()
            .unwrap_or_else(|error| panic!("failed to start {:?}: {error}", python));
        assert!(
            status.success(),
            "failed to compile embedded Python asset {}",
            source.display()
        );
    }
}
