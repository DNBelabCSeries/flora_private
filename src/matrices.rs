use std::collections::BTreeSet;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

pub fn is_genomic_placeholder(gene: &str) -> bool {
    let mut parts = gene.rsplitn(3, '_');
    let Some(c) = parts.next() else { return false };
    let Some(b) = parts.next() else { return false };
    let Some(a) = parts.next() else { return false };
    !a.is_empty()
        && !b.is_empty()
        && !c.is_empty()
        && b.chars().all(|ch| ch.is_ascii_digit())
        && c.chars().all(|ch| ch.is_ascii_digit())
}

pub fn add_unique_umi(
    matrix: &mut HashMap<(String, String), HashSet<String>>,
    feature: &str,
    cell: &str,
    umi: &str,
) {
    matrix
        .entry((feature.to_string(), cell.to_string()))
        .or_default()
        .insert(umi.to_string());
}

pub fn matrix_axes(
    matrix: &HashMap<(String, String), HashSet<String>>,
) -> (BTreeSet<String>, BTreeSet<String>) {
    let mut rows = BTreeSet::new();
    let mut cols = BTreeSet::new();
    for (feature, cell) in matrix.keys() {
        rows.insert(feature.clone());
        cols.insert(cell.clone());
    }
    (rows, cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genomic_placeholder_detection_matches_python_intent() {
        assert!(is_genomic_placeholder("chr7_44468000_44469000"));
        assert!(is_genomic_placeholder("hs_chr7_44468000_44469000"));
        assert!(is_genomic_placeholder("NC_000001.11_1000_2000"));
        assert!(!is_genomic_placeholder("ACTB"));
        assert!(!is_genomic_placeholder("chr7_44468000_gene"));
    }

    #[test]
    fn unique_umi_counts_deduplicate_per_feature_cell() {
        let mut matrix = HashMap::default();
        add_unique_umi(&mut matrix, "gene1", "cell1", "umi1");
        add_unique_umi(&mut matrix, "gene1", "cell1", "umi1");
        add_unique_umi(&mut matrix, "gene1", "cell1", "umi2");
        assert_eq!(matrix[&("gene1".to_string(), "cell1".to_string())].len(), 2);
    }
}
