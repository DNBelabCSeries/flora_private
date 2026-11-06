use std::collections::VecDeque;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use crate::barcode::bounded_levenshtein;

pub fn cluster_directional(
    counts: &HashMap<String, usize>,
    threshold: usize,
) -> HashMap<String, String> {
    let graph = adjacency(counts, threshold);
    let components = connected_components(&graph, counts);
    create_map_to_correct_umi(components)
}

fn adjacency(counts: &HashMap<String, usize>, threshold: usize) -> HashMap<String, Vec<String>> {
    let mut adj = counts
        .keys()
        .cloned()
        .map(|umi| (umi, Vec::new()))
        .collect::<HashMap<_, _>>();
    let umis = counts.keys().cloned().collect::<Vec<_>>();
    for (i, j) in candidate_pairs(&umis, threshold) {
        if i >= j {
            continue;
        }
        let a = &umis[i];
        let b = &umis[j];
        if bounded_levenshtein(a, b, threshold).is_none() {
            continue;
        }
        let ca = counts[a];
        let cb = counts[b];
        if ca >= (cb * 2).saturating_sub(1) {
            if let Some(neighbors) = adj.get_mut(a) {
                neighbors.push(b.clone());
            }
        }
        if cb >= (ca * 2).saturating_sub(1) {
            if let Some(neighbors) = adj.get_mut(b) {
                neighbors.push(a.clone());
            }
        }
    }
    adj
}

fn candidate_pairs(umis: &[String], threshold: usize) -> Vec<(usize, usize)> {
    if umis.len() < 2 {
        return Vec::new();
    }
    if threshold == 0 || umis.len() <= 64 {
        return all_pairs(umis.len());
    }
    let same_len = umis
        .first()
        .map(|umi| umi.len())
        .map(|len| umis.iter().all(|umi| umi.len() == len))
        .unwrap_or(true);
    if !same_len {
        return all_pairs(umis.len());
    }

    let len = umis[0].len();
    let segments = threshold + 1;
    if len < segments {
        return all_pairs(umis.len());
    }

    let bounds = segment_bounds(len, segments);
    let mut buckets: HashMap<(usize, Vec<u8>), Vec<usize>> = HashMap::default();
    for (umi_idx, umi) in umis.iter().enumerate() {
        let bytes = umi.as_bytes();
        for (segment_idx, (start, end)) in bounds.iter().copied().enumerate() {
            buckets
                .entry((segment_idx, bytes[start..end].to_vec()))
                .or_default()
                .push(umi_idx);
        }
    }

    let mut seen_pairs = HashSet::default();
    let mut pairs = Vec::new();
    // One segment survives threshold edits intact, but indels can shift it.
    for (i, umi) in umis.iter().enumerate() {
        let bytes = umi.as_bytes();
        for (segment_idx, (start, end)) in bounds.iter().copied().enumerate() {
            let segment_len = end - start;
            let first = start.saturating_sub(threshold);
            let last = start.saturating_add(threshold).min(len - segment_len);
            for offset in first..=last {
                let key = (segment_idx, bytes[offset..offset + segment_len].to_vec());
                let Some(indices) = buckets.get(&key) else {
                    continue;
                };
                for &j in indices {
                    if i == j {
                        continue;
                    }
                    let pair = if i < j { (i, j) } else { (j, i) };
                    let packed = ((pair.0 as u64) << 32) | pair.1 as u64;
                    if seen_pairs.insert(packed) {
                        pairs.push(pair);
                    }
                }
            }
        }
    }
    pairs
}

fn all_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
    for i in 0..n {
        for j in i + 1..n {
            pairs.push((i, j));
        }
    }
    pairs
}

fn segment_bounds(len: usize, segments: usize) -> Vec<(usize, usize)> {
    let mut bounds = Vec::with_capacity(segments);
    for idx in 0..segments {
        let start = idx * len / segments;
        let end = (idx + 1) * len / segments;
        bounds.push((start, end));
    }
    bounds
}

fn connected_components(
    graph: &HashMap<String, Vec<String>>,
    counts: &HashMap<String, usize>,
) -> Vec<Vec<String>> {
    let mut nodes = graph.keys().cloned().collect::<Vec<_>>();
    nodes.sort_by(|a, b| counts[b].cmp(&counts[a]).then_with(|| a.cmp(b)));
    let mut seen = HashSet::default();
    let mut components = Vec::new();

    for node in nodes {
        if !seen.insert(node.clone()) {
            continue;
        }
        let mut queue = VecDeque::from([node.clone()]);
        let mut component = vec![node];
        while let Some(curr) = queue.pop_front() {
            if let Some(neighbors) = graph.get(&curr) {
                for next in neighbors {
                    if seen.insert(next.clone()) {
                        queue.push_back(next.clone());
                        component.push(next.clone());
                    }
                }
            }
        }
        component.sort_by(|a, b| counts[b].cmp(&counts[a]).then_with(|| a.cmp(b)));
        components.push(component);
    }
    components
}

fn create_map_to_correct_umi(components: Vec<Vec<String>>) -> HashMap<String, String> {
    let mut map = HashMap::default();
    for component in components {
        if let Some(head) = component.first().cloned() {
            for umi in component {
                map.insert(umi, head.clone());
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directional_clustering_merges_low_count_neighbors() {
        let counts = [
            ("AAAA".to_string(), 10usize),
            ("AAAT".to_string(), 3usize),
            ("TTTT".to_string(), 2usize),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        let map = cluster_directional(&counts, 1);
        assert_eq!(map["AAAA"], "AAAA");
        assert_eq!(map["AAAT"], "AAAA");
        assert_eq!(map["TTTT"], "TTTT");
    }

    #[test]
    fn candidate_pairs_keeps_edit_distance_neighbors_via_segment_buckets() {
        let umis = vec![
            "ACGTACGT".to_string(),
            "ACGTTCGT".to_string(),
            "TTTTTTTT".to_string(),
        ];
        let pairs = candidate_pairs(&umis, 1);
        assert!(pairs.contains(&(0, 1)));
    }

    #[test]
    fn indexed_candidates_keep_shifted_indel_neighbors() {
        let mut umis = vec![
            "ACGTCAGTACGATCGTACGA".to_string(),
            "CGTCAGTACGATCGTACGAA".to_string(),
        ];
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            umis.push(format!("TTTTTTTTTTTT{suffix}"));
        }
        assert_eq!(bounded_levenshtein(&umis[0], &umis[1], 2), Some(2));
        for threshold in 2..=3 {
            assert!(candidate_pairs(&umis, threshold).contains(&(0, 1)));
        }
    }

    #[test]
    fn clustering_keeps_indel_correction_across_indexing_boundary() {
        let parent = "ACGTCAGTACGATCGTACGA";
        let child = "CGTCAGTACGATCGTACGAA";
        for distinct_umis in [64, 65] {
            let mut counts = HashMap::default();
            counts.insert(parent.to_string(), 100);
            counts.insert(child.to_string(), 1);
            for i in 0..distinct_umis - 2 {
                let suffix = (0..8)
                    .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                    .collect::<String>();
                counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
            }
            assert_eq!(cluster_directional(&counts, 2)[child], parent);
        }
    }

    #[test]
    fn indexed_clustering_rejects_distant_segment_matches() {
        let parent = "AAAAAAAAAAAAAAAAAAAA";
        let distant = "AAAAAAAACCCCCCCCCCCC";
        let mut counts = HashMap::default();
        counts.insert(parent.to_string(), 100);
        counts.insert(distant.to_string(), 1);
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
        }
        let map = cluster_directional(&counts, 2);
        assert_eq!(map[parent], parent);
        assert_eq!(map[distant], distant);
    }

    #[test]
    fn indexed_clustering_preserves_directional_abundance_filter() {
        let parent = "ACGTCAGTACGATCGTACGA";
        let child = "CGTCAGTACGATCGTACGAA";
        let mut counts = HashMap::default();
        counts.insert(parent.to_string(), 100);
        counts.insert(child.to_string(), 51);
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
        }
        let map = cluster_directional(&counts, 2);
        assert_eq!(map[parent], parent);
        assert_eq!(map[child], child);
    }

    #[test]
    fn indexed_candidates_match_exhaustive_binary_dna_neighbors() {
        let umis = (0..256)
            .map(|value| {
                (0..8)
                    .map(|bit| if value & (1 << bit) == 0 { 'A' } else { 'C' })
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        for threshold in 1..=3 {
            let candidates = candidate_pairs(&umis, threshold)
                .into_iter()
                .collect::<HashSet<_>>();
            for (i, j) in all_pairs(umis.len()) {
                if bounded_levenshtein(&umis[i], &umis[j], threshold).is_some() {
                    assert!(
                        candidates.contains(&(i, j)),
                        "missing {i}/{j}, distance <= {threshold}"
                    );
                }
            }
        }
    }
}
