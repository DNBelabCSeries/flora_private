use std::cmp::min;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use serde::Serialize;

use crate::fastq::FastqRecord;

#[derive(Clone, Debug, Default, Serialize)]
pub struct PutativeRow {
    pub read_id: String,
    pub putative_bc: String,
    pub bc_fixed_locs: isize,
    pub putative_bc_min_qs: Option<i32>,
    pub putative_umi: String,
    pub umi_fixed_locs: Option<isize>,
    pub post_umi_flankings: String,
    pub poly_a_starts: Option<isize>,
    pub read_types: u8,
    pub putative_bc_5p: String,
    pub bc_fixed_locs_50: isize,
    pub putative_bc_min_qs_5p: Option<i32>,
    pub putative_umi_5p: String,
    pub umi_fixed_locs_5p: Option<isize>,
}

#[derive(Clone, Debug)]
pub struct CorrectedRead {
    pub read_id: String,
    pub putative_umi: String,
    pub putative_umi_5p: String,
    pub bc3_corrected: String,
    pub bc5_corrected: String,
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            other => other,
        })
        .collect()
}

pub fn revcomp_upper(seq: &str) -> String {
    seq.trim()
        .to_ascii_uppercase()
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

fn slice(seq: &str, start: isize, end: isize) -> String {
    let n = seq.len() as isize;
    let s = if start < 0 { n + start } else { start }.clamp(0, n) as usize;
    let e = if end < 0 { n + end } else { end }.clamp(0, n) as usize;
    if s >= e {
        String::new()
    } else {
        seq[s..e].to_string()
    }
}

fn min_q(qual: &str) -> Option<i32> {
    qual.bytes().map(|b| b as i32 - 33).min()
}

fn rfind_with_negative(s: &str, sub: &str) -> isize {
    s.rfind(sub)
        .map(|pos| pos as isize - s.len() as isize)
        .unwrap_or(-1)
}

fn find_pos(s: &str, sub: &str) -> isize {
    s.find(sub).map(|pos| pos as isize).unwrap_or(-1)
}

pub fn extract_putative(
    rec: &FastqRecord,
    bc_fixed_3p: &str,
    umi_fixed_3p: &str,
    bc_fixed_5p: &str,
    umi_fixed_5p: &str,
) -> PutativeRow {
    let (bc3, bc3_loc, q3, umi3, umi3_loc, flanking3, poly_a_start, read_type) =
        extract_3p(rec, bc_fixed_3p, umi_fixed_3p);
    let (bc5, bc5_loc, q5, umi5, umi5_loc) = extract_5p(rec, bc_fixed_5p, umi_fixed_5p);
    PutativeRow {
        read_id: rec.id.clone(),
        putative_bc: bc3,
        bc_fixed_locs: bc3_loc,
        putative_bc_min_qs: q3,
        putative_umi: umi3,
        umi_fixed_locs: umi3_loc,
        post_umi_flankings: flanking3,
        poly_a_starts: poly_a_start,
        read_types: read_type,
        putative_bc_5p: bc5,
        bc_fixed_locs_50: bc5_loc,
        putative_bc_min_qs_5p: q5,
        putative_umi_5p: umi5,
        umi_fixed_locs_5p: umi5_loc,
    }
}

fn absolute_index(index: isize, read_len: usize) -> Option<usize> {
    let index = if index < 0 {
        read_len as isize + index
    } else {
        index
    };
    (0..=read_len as isize)
        .contains(&index)
        .then_some(index as usize)
}

/// Return the cDNA interval after removing recognized terminal barcode/UMI
/// structures. The booleans come from the corrected read, so an unvalidated
/// end is never trimmed solely because a short motif happened to match.
pub fn cdna_trim_bounds(
    row: &PutativeRow,
    read_len: usize,
    has_valid_3p: bool,
    has_valid_5p: bool,
) -> Option<(usize, usize)> {
    let mut start = 0usize;
    let mut end = read_len;

    if has_valid_5p {
        if let Some(umi_fixed) = row
            .umi_fixed_locs_5p
            .and_then(|index| absolute_index(index, read_len))
        {
            start = umi_fixed.saturating_add(15).min(read_len);
        } else if row.bc_fixed_locs_50 >= 0 {
            start = (row.bc_fixed_locs_50 as usize)
                .saturating_add(31)
                .min(read_len);
        }
    }

    if has_valid_3p {
        if let Some(poly_a_start) = row
            .poly_a_starts
            .and_then(|index| absolute_index(index, read_len))
        {
            end = poly_a_start;
        } else if let Some(umi_fixed) = row
            .umi_fixed_locs
            .and_then(|index| absolute_index(index, read_len))
        {
            end = umi_fixed.saturating_sub(10);
        } else if row.bc_fixed_locs != -1 {
            if let Some(bc_fixed) = absolute_index(row.bc_fixed_locs, read_len) {
                end = bc_fixed.saturating_sub(25);
            }
        }
    }

    (start < end && end <= read_len).then_some((start, end))
}

fn poly_a_trimming_idx(
    seq: &str,
    seed: &str,
    window: usize,
    min_a: usize,
    min_tail_len: usize,
) -> Option<usize> {
    let s = seq.to_ascii_uppercase();
    let anchor = s.rfind(seed)?;
    let mut poly_a_start = anchor;
    let mut i = anchor as isize - 1;
    while i >= 0 {
        let idx = i as usize;
        if s.as_bytes()[idx] == b'A' {
            poly_a_start = idx;
            i -= 1;
            continue;
        }
        let left = idx.saturating_sub(window.saturating_sub(1));
        let a_count = s.as_bytes()[left..=idx]
            .iter()
            .filter(|b| **b == b'A')
            .count();
        if a_count >= min_a {
            poly_a_start = idx;
            i -= 1;
            continue;
        }
        break;
    }
    (s.len().saturating_sub(poly_a_start) >= min_tail_len).then_some(poly_a_start)
}

fn poly_a_trimming_idx_neg(seq: &str) -> Option<isize> {
    poly_a_trimming_idx(seq, "AAAA", 10, 7, 8).map(|idx| idx as isize - seq.len() as isize)
}

fn poly_a_start_before(rec: &FastqRecord, end: isize) -> Option<isize> {
    let seq_poly_a = slice(&rec.seq, end - 100, end);
    poly_a_trimming_idx_neg(&seq_poly_a).map(|idx| idx + end)
}

fn extract_3p(
    rec: &FastqRecord,
    bc_fixed: &str,
    umi_fixed: &str,
) -> (
    String,
    isize,
    Option<i32>,
    String,
    Option<isize>,
    String,
    Option<isize>,
    u8,
) {
    let part_seq = slice(&rec.seq, -30, rec.seq.len() as isize);
    let bc_loc = rfind_with_negative(&part_seq, bc_fixed);
    let bc: String;
    let mut q = None;
    let mut umi = String::new();
    let mut umi_loc = None;
    let mut flank = String::new();
    let mut poly_a_start = None;
    let read_type;

    if bc_loc == -16 {
        bc = slice(&rec.seq, -26, rec.seq.len() as isize);
        q = min_q(&slice(&rec.qual, -26, rec.qual.len() as isize));
        let find_umi_seq = slice(&rec.seq, -36, -26);
        let urel = rfind_with_negative(&find_umi_seq, umi_fixed);
        if urel != -1 {
            let loc = urel - 26;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc - 10, loc + 5);
            flank = slice(&rec.seq, loc - 15, loc - 10);
            poly_a_start = poly_a_start_before(rec, loc - 10);
            read_type = if poly_a_start.is_some() { 1 } else { 2 };
        } else {
            umi = slice(&rec.seq, -41, -26);
            flank = slice(&rec.seq, -46, -41);
            poly_a_start = poly_a_start_before(rec, -26);
            read_type = if poly_a_start.is_some() { 3 } else { 4 };
        }
    } else if bc_loc < -16 {
        bc = slice(&rec.seq, bc_loc - 10, bc_loc + 16);
        q = min_q(&slice(&rec.qual, bc_loc - 10, bc_loc + 16));
        let find_umi_seq = slice(&rec.seq, bc_loc - 20, bc_loc - 10);
        let urel = rfind_with_negative(&find_umi_seq, umi_fixed);
        if urel != -1 {
            let loc = urel + bc_loc - 10;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc - 10, loc + 5);
            flank = slice(&rec.seq, loc - 15, loc - 10);
            poly_a_start = poly_a_start_before(rec, loc - 10);
            read_type = if poly_a_start.is_some() { 5 } else { 6 };
        } else {
            umi = slice(&rec.seq, bc_loc - 25, bc_loc - 10);
            flank = slice(&rec.seq, bc_loc - 30, bc_loc - 25);
            poly_a_start = poly_a_start_before(rec, bc_loc - 10);
            read_type = if poly_a_start.is_some() { 7 } else { 8 };
        }
    } else if bc_loc == -1 {
        bc = slice(&rec.seq, -26, rec.seq.len() as isize);
        let find_umi_seq = slice(&rec.seq, 0, -40);
        let urel = rfind_with_negative(&find_umi_seq, umi_fixed);
        if urel != -1 {
            umi_loc = Some(urel);
            umi = slice(&rec.seq, urel - 10, urel + 5);
            flank = slice(&rec.seq, urel - 15, urel - 10);
            poly_a_start = poly_a_start_before(rec, urel - 10);
            read_type = if poly_a_start.is_some() { 9 } else { 10 };
        } else {
            read_type = 11;
        }
    } else {
        let find_umi_seq = slice(&rec.seq, bc_loc - 20, bc_loc - 10);
        let urel = rfind_with_negative(&find_umi_seq, umi_fixed);
        if urel != -1 {
            let loc = urel + bc_loc - 10;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc - 10, loc + 5);
            flank = slice(&rec.seq, loc - 15, loc - 10);
            let start = loc + 5;
            bc = slice(&rec.seq, start, rec.seq.len() as isize);
            q = min_q(&slice(&rec.qual, start, rec.qual.len() as isize));
            poly_a_start = poly_a_start_before(rec, loc - 10);
            read_type = if poly_a_start.is_some() { 12 } else { 13 };
        } else {
            bc = slice(&rec.seq, -26, rec.seq.len() as isize);
            read_type = 14;
        }
    }
    (bc, bc_loc, q, umi, umi_loc, flank, poly_a_start, read_type)
}

fn extract_5p(
    rec: &FastqRecord,
    bc_fixed: &str,
    umi_fixed: &str,
) -> (String, isize, Option<i32>, String, Option<isize>) {
    let part_seq = slice(&rec.seq, 0, 30);
    let bc_loc = find_pos(&part_seq, bc_fixed);
    let bc;
    let q;
    let mut umi = String::new();
    let mut umi_loc = None;

    if bc_loc == 11 {
        bc = slice(&rec.seq, 0, 26);
        q = min_q(&slice(&rec.qual, 0, 26));
        let find_umi_seq = slice(&rec.seq, 26, 36);
        let urel = find_pos(&find_umi_seq, umi_fixed);
        if urel != -1 {
            let loc = urel + 26;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc, loc + 15);
        } else {
            umi = slice(&rec.seq, 26, 41);
        }
    } else if bc_loc > 11 {
        bc = slice(&rec.seq, bc_loc - 10, bc_loc + 16);
        q = min_q(&slice(&rec.qual, bc_loc - 10, bc_loc + 16));
        let start = bc_loc + 16;
        let urel = find_pos(&slice(&rec.seq, start, start + 10), umi_fixed);
        if urel != -1 {
            let loc = urel + start;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc, loc + 15);
        } else {
            umi = slice(&rec.seq, start, start + 15);
        }
    } else if bc_loc >= 0 {
        bc = slice(&rec.seq, 0, bc_loc + 16);
        q = min_q(&slice(&rec.qual, 0, bc_loc + 16));
        let start = bc_loc + 16;
        let urel = find_pos(&slice(&rec.seq, start, start + 10), umi_fixed);
        if urel != -1 {
            let loc = urel + start;
            umi_loc = Some(loc);
            umi = slice(&rec.seq, loc, loc + 15);
        } else {
            umi = slice(&rec.seq, start, start + 15);
        }
    } else {
        bc = slice(&rec.seq, 0, 26);
        q = None;
        let urel = find_pos(&slice(&rec.seq, 0, 40), umi_fixed);
        if urel != -1 {
            umi_loc = Some(urel);
            umi = slice(&rec.seq, urel, urel + 15);
        }
    }
    (bc, bc_loc, q, umi, umi_loc)
}

pub fn bounded_levenshtein(a: &str, b: &str, max_ed: usize) -> Option<usize> {
    if a.len().abs_diff(b.len()) > max_ed {
        return None;
    }
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0; b.len() + 1];
    for (i, ca) in a.bytes().enumerate() {
        curr[0] = i + 1;
        let mut row_min = curr[0];
        for (j, cb) in b.bytes().enumerate() {
            let cost = usize::from(ca != cb);
            curr[j + 1] = min(min(curr[j] + 1, prev[j + 1] + 1), prev[j] + cost);
            row_min = row_min.min(curr[j + 1]);
        }
        if row_min > max_ed {
            return None;
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    (prev[b.len()] <= max_ed).then_some(prev[b.len()])
}

fn bounded_sub_edit_distance(pattern: &str, text: &str, max_ed: usize) -> Option<(usize, usize)> {
    if pattern.is_empty() {
        return Some((0, 0));
    }
    let m = pattern.len();
    let n = text.len();
    let p = pattern.as_bytes();
    let t = text.as_bytes();

    let mut prev = vec![0usize; n + 1];
    let mut curr = vec![0usize; n + 1];
    for i in 1..=m {
        curr[0] = i;
        let mut row_min = curr[0];
        for j in 1..=n {
            let cost = usize::from(p[i - 1] != t[j - 1]);
            curr[j] = min(min(curr[j - 1] + 1, prev[j] + 1), prev[j - 1] + cost);
            row_min = row_min.min(curr[j]);
        }
        if row_min > max_ed {
            return None;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let mut best = usize::MAX;
    let mut best_end = 0usize;
    for (end, ed) in prev.iter().copied().enumerate() {
        if ed < best {
            best = ed;
            best_end = end;
        }
    }
    (best <= max_ed).then_some((best, best_end))
}

#[derive(Clone, Debug, Default)]
pub struct BarcodeIndex {
    exact: HashSet<String>,
    values: Vec<String>,
    segment_index: HashMap<String, Vec<usize>>,
    segment_lengths: Vec<usize>,
    indexed_max_ed: Option<usize>,
}

impl BarcodeIndex {
    pub fn new<I>(values: I) -> Self
    where
        I: IntoIterator<Item = String>,
    {
        Self::from_values(values, None)
    }

    pub fn with_max_ed<I>(values: I, max_ed: usize) -> Self
    where
        I: IntoIterator<Item = String>,
    {
        Self::from_values(values, Some(max_ed))
    }

    pub fn contains(&self, value: &str) -> bool {
        self.exact.contains(value)
    }

    fn from_values<I>(values: I, indexed_max_ed: Option<usize>) -> Self
    where
        I: IntoIterator<Item = String>,
    {
        let mut index = Self {
            indexed_max_ed,
            ..Self::default()
        };
        for value in values {
            if !value.is_empty() && index.exact.insert(value.clone()) {
                index.values.push(value);
            }
        }
        if let Some(max_ed) = indexed_max_ed {
            index.build_segment_index(max_ed);
        }
        index
    }

    fn build_segment_index(&mut self, max_ed: usize) {
        let parts = max_ed.saturating_add(1);
        if parts == 0 {
            return;
        }
        let mut lengths = HashSet::default();
        for (candidate_id, candidate) in self.values.iter().enumerate() {
            if parts > candidate.len() {
                self.segment_index.clear();
                self.segment_lengths.clear();
                self.indexed_max_ed = None;
                return;
            }
            for (start, end) in partition_ranges(candidate.len(), parts) {
                let segment = candidate[start..end].to_string();
                lengths.insert(end - start);
                self.segment_index
                    .entry(segment)
                    .or_default()
                    .push(candidate_id);
            }
        }
        self.segment_lengths = lengths.into_iter().collect();
        self.segment_lengths.sort_unstable();
    }

    fn candidate_ids(&self, query: &str, max_ed: usize) -> Vec<usize> {
        if self.indexed_max_ed != Some(max_ed) || self.segment_index.is_empty() {
            return (0..self.values.len()).collect();
        }
        let mut hits = HashSet::default();
        for &length in &self.segment_lengths {
            if length > query.len() {
                continue;
            }
            for start in 0..=query.len() - length {
                if let Some(ids) = self.segment_index.get(&query[start..start + length]) {
                    hits.extend(ids.iter().copied());
                }
            }
        }
        let mut hits = hits.into_iter().collect::<Vec<_>>();
        hits.sort_unstable();
        hits
    }
}

fn partition_ranges(length: usize, parts: usize) -> Vec<(usize, usize)> {
    let base = length / parts;
    let remainder = length % parts;
    let mut start = 0usize;
    (0..parts)
        .map(|part| {
            let width = base + usize::from(part < remainder);
            let range = (start, start + width);
            start += width;
            range
        })
        .collect()
}

pub fn correct_one_side(bc: &str, whitelist: &HashSet<String>, max_ed: usize) -> String {
    if bc.trim().is_empty() {
        return String::new();
    }
    if whitelist.contains(bc) {
        return bc.to_string();
    }
    let mut best_ed = max_ed;
    let mut bc_hit: Option<&String> = None;
    let mut ambiguous = false;
    for candidate in whitelist {
        if let Some((ed, _)) = bounded_sub_edit_distance(candidate, bc, best_ed) {
            if ed < best_ed {
                best_ed = ed;
                bc_hit = Some(candidate);
                ambiguous = false;
            } else if ed == best_ed {
                if bc_hit.is_none() {
                    bc_hit = Some(candidate);
                    ambiguous = false;
                } else {
                    ambiguous = true;
                    if best_ed == 0 {
                        return String::new();
                    }
                    best_ed -= 1;
                }
            }
        }
    }
    if ambiguous {
        String::new()
    } else {
        bc_hit.cloned().unwrap_or_default()
    }
}

pub fn correct_one_side_indexed(bc: &str, index: &BarcodeIndex, max_ed: usize) -> String {
    if bc.trim().is_empty() {
        return String::new();
    }
    if index.contains(bc) {
        return bc.to_string();
    }
    let mut best_ed = max_ed;
    let mut bc_hit: Option<&String> = None;
    let mut ambiguous = false;
    for candidate_id in index.candidate_ids(bc, max_ed) {
        let candidate = &index.values[candidate_id];
        if let Some((ed, _)) = bounded_sub_edit_distance(candidate, bc, best_ed) {
            if ed < best_ed {
                best_ed = ed;
                bc_hit = Some(candidate);
                ambiguous = false;
            } else if ed == best_ed {
                if bc_hit.is_none() {
                    bc_hit = Some(candidate);
                    ambiguous = false;
                } else {
                    ambiguous = true;
                    if best_ed == 0 {
                        return String::new();
                    }
                    best_ed -= 1;
                }
            }
        }
    }
    if ambiguous {
        String::new()
    } else {
        bc_hit.cloned().unwrap_or_default()
    }
}

pub fn correct_dual(
    row: &PutativeRow,
    wl3: &HashSet<String>,
    wl5: &HashSet<String>,
    max_ed: usize,
    min_q: i32,
) -> CorrectedRead {
    let mut bc3 = row.putative_bc.clone();
    let mut bc5 = row.putative_bc_5p.clone();
    let mut umi3 = row.putative_umi.clone();
    let mut umi5 = row.putative_umi_5p.clone();
    if row.putative_bc_min_qs.is_some_and(|q| q < min_q) {
        bc3.clear();
        umi3.clear();
    }
    if row.putative_bc_min_qs_5p.is_some_and(|q| q < min_q) {
        bc5.clear();
        umi5.clear();
    }
    let c3 = if wl3.contains(&bc3) {
        bc3
    } else {
        correct_one_side(&bc3, wl3, max_ed)
    };
    let c5 = if wl5.contains(&bc5) {
        bc5
    } else {
        correct_one_side(&bc5, wl5, max_ed)
    };
    CorrectedRead {
        read_id: row.read_id.clone(),
        putative_umi: if c3.is_empty() { String::new() } else { umi3 },
        putative_umi_5p: if c5.is_empty() { String::new() } else { umi5 },
        bc3_corrected: c3,
        bc5_corrected: c5,
    }
}

fn strip_fixed(seq: &str, middle6: &str) -> String {
    let s = seq.trim().to_ascii_uppercase();
    if s.len() == 26 && &s[10..16] == middle6 {
        format!("{}{}", &s[..10], &s[16..])
    } else {
        String::new()
    }
}

pub fn strip_fixed_3p(seq: &str) -> String {
    strip_fixed(seq, "GCTACC")
}

pub fn strip_fixed_5p(seq: &str) -> String {
    strip_fixed(seq, "CCTTCC")
}

pub fn umi_a_ratio(umi: &str) -> f64 {
    let s = umi.trim().to_ascii_uppercase();
    if s.is_empty() {
        return 0.0;
    }
    s.bytes().filter(|b| *b == b'A').count() as f64 / s.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_matches_python_orientation() {
        assert_eq!(reverse_complement("GGTAGC"), "GCTACC");
        assert_eq!(revcomp_upper("ACGTN"), "NACGT");
    }

    #[test]
    fn bounded_edit_distance_respects_limit() {
        assert_eq!(bounded_levenshtein("AAAA", "AAAT", 1), Some(1));
        assert_eq!(bounded_levenshtein("AAAA", "TTTT", 2), None);
    }

    #[test]
    fn correction_rejects_ambiguous_best_hit() {
        let wl = ["AAAA".to_string(), "AAAT".to_string()]
            .into_iter()
            .collect::<HashSet<_>>();
        assert_eq!(correct_one_side("AAAC", &wl, 2), "");
    }

    #[test]
    fn indexed_correction_finds_nearby_barcode() {
        let index = BarcodeIndex::new(["AAAACCCC".to_string(), "TTTTGGGG".to_string()]);
        assert_eq!(correct_one_side_indexed("AAAACCCA", &index, 1), "AAAACCCC");
        assert_eq!(
            correct_one_side_indexed("TAAAACCCCA", &index, 1),
            "AAAACCCC"
        );
        assert_eq!(correct_one_side_indexed("CCCCCCCC", &index, 1), "");
    }

    #[test]
    fn segmented_index_matches_exhaustive_semiglobal_correction() {
        let mut seed = 0x5eed_u64;
        let mut whitelist = HashSet::default();
        for _ in 0..250 {
            whitelist.insert(random_dna(&mut seed, 26));
        }
        let index = BarcodeIndex::with_max_ed(whitelist.iter().cloned(), 3);
        for _ in 0..500 {
            let length = 22 + (next_u64(&mut seed) % 10) as usize;
            let query = random_dna(&mut seed, length);
            for max_ed in 1..=3 {
                assert_eq!(
                    correct_one_side_indexed(&query, &index, max_ed),
                    correct_one_side(&query, &whitelist, max_ed),
                    "query={query} max_ed={max_ed}"
                );
            }
        }
    }

    #[test]
    fn segmented_index_preserves_free_flanks_indels_and_ambiguity() {
        let whitelist = [
            "AAAACCCCGGGGTTTT".to_string(),
            "AAAACCCCGGGGTTTA".to_string(),
            "TTTTGGGGCCCCAAAA".to_string(),
        ]
        .into_iter()
        .collect::<HashSet<_>>();
        let index = BarcodeIndex::with_max_ed(whitelist.iter().cloned(), 2);
        for query in [
            "GGAAAACCCCGGGGTTTTCC",
            "AAAACCCGGGGTTTT",
            "AAAACCCCCGGGGTTTT",
            "AAAACCCCGGGGTTTC",
        ] {
            assert_eq!(
                correct_one_side_indexed(query, &index, 2),
                correct_one_side(query, &whitelist, 2),
                "query={query}"
            );
        }
    }

    #[test]
    fn segmented_index_exhaustively_matches_small_binary_alphabet() {
        let whitelist = all_binary_dna(4).into_iter().collect::<HashSet<_>>();
        let index = BarcodeIndex::with_max_ed(whitelist.iter().cloned(), 2);
        for length in 1..=7 {
            for query in all_binary_dna(length) {
                for max_ed in 1..=2 {
                    assert_eq!(
                        correct_one_side_indexed(&query, &index, max_ed),
                        correct_one_side(&query, &whitelist, max_ed),
                        "query={query} max_ed={max_ed}"
                    );
                }
            }
        }
    }

    fn next_u64(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        *seed
    }

    fn random_dna(seed: &mut u64, length: usize) -> String {
        (0..length)
            .map(|_| match next_u64(seed) & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            })
            .collect()
    }

    fn all_binary_dna(length: usize) -> Vec<String> {
        (0..(1usize << length))
            .map(|bits| {
                (0..length)
                    .map(|shift| if (bits >> shift) & 1 == 0 { 'A' } else { 'C' })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn extracts_exact_5p_and_3p_barcodes() {
        let seq =
            "AAAAAAAAAACCTTCCGGGGGGGGGGCAGCATTTTTTTTTTAAAAAAAAGATCTCCCCCCCCCCGCTACCTTTTTTTTTT";
        let rec = FastqRecord {
            id: "read1".to_string(),
            seq: seq.to_string(),
            qual: "I".repeat(seq.len()),
        };
        let row = extract_putative(&rec, "GCTACC", "AGATC", "CCTTCC", "CAGCA");
        assert_eq!(row.putative_bc_5p, "AAAAAAAAAACCTTCCGGGGGGGGGG");
        assert_eq!(row.putative_bc, "CCCCCCCCCCGCTACCTTTTTTTTTT");
        assert_eq!(row.putative_umi_5p, "CAGCATTTTTTTTTT");
        assert_eq!(row.putative_umi, "TTTAAAAAAAAGATC");
        assert_eq!(row.read_types, 1);
        assert!(row.poly_a_starts.is_some());
    }

    #[test]
    fn trims_valid_dual_end_structures_to_cdna() {
        let barcode_5p = "AAAAAAAAAACCTTCCGGGGGGGGGG";
        let umi_5p = "CAGCATTTTTTTTTT";
        let cdna = "CGTC".repeat(40);
        let umi_3p = "CCCCCCCCCCAGATC";
        let barcode_3p = "TTTTTTTTTTGCTACCAAAAAAAAAA";
        let seq = format!("{barcode_5p}{umi_5p}{cdna}{umi_3p}{barcode_3p}");
        let rec = FastqRecord {
            id: "trim-me".to_string(),
            qual: "I".repeat(seq.len()),
            seq,
        };
        let row = extract_putative(&rec, "GCTACC", "AGATC", "CCTTCC", "CAGCA");
        let (start, end) = cdna_trim_bounds(&row, rec.seq.len(), true, true).unwrap();
        assert_eq!(&rec.seq[start..end], cdna);
        assert_eq!(end - start, cdna.len());
    }

    #[test]
    fn leaves_unvalidated_ends_untouched() {
        let rec = FastqRecord {
            id: "no-trim".to_string(),
            seq: "ACGTACGT".to_string(),
            qual: "IIIIIIII".to_string(),
        };
        let row = extract_putative(&rec, "GCTACC", "AGATC", "CCTTCC", "CAGCA");
        assert_eq!(
            cdna_trim_bounds(&row, rec.seq.len(), false, false),
            Some((0, 8))
        );
    }
}
