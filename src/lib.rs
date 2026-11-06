extern crate self as flora;

pub mod annotation;
pub mod bam_runtime;
pub mod barcode;
pub mod fastq;
pub mod glycine;
pub mod matrices;
pub mod pipeline;
pub mod read_qc;
mod trim_cache;
pub mod umi_cluster;
pub mod workflow_runtime;

#[path = "bin/add_cb_ur_tags.rs"]
pub mod stage_add_cb_ur_tags;
#[path = "bin/add_gene_tags.rs"]
pub mod stage_add_gene_tags;
#[path = "bin/assign_genes.rs"]
pub mod stage_assign_genes;
#[path = "bin/assign_transcripts.rs"]
pub mod stage_assign_transcripts;
#[path = "bin/cell_umi_gene_table.rs"]
pub mod stage_cell_umi_gene_table;
#[path = "bin/cluster_umis_allbam.rs"]
pub mod stage_cluster_umis_allbam;
#[path = "bin/gene_expression.rs"]
pub mod stage_gene_expression;
#[path = "bin/generate_26bp_whitelists.rs"]
pub mod stage_generate_26bp_whitelists;
#[path = "bin/isoform_expression.rs"]
pub mod stage_isoform_expression;
#[path = "bin/prepare_read_tags.rs"]
pub mod stage_prepare_read_tags;
#[path = "bin/read_qc_summary.rs"]
pub mod stage_read_qc_summary;
#[path = "bin/rna_qc_metrics.rs"]
pub mod stage_rna_qc_metrics;
#[path = "bin/saturation.rs"]
pub mod stage_saturation;
#[path = "bin/tag_and_assign_genes.rs"]
pub mod stage_tag_and_assign_genes;
