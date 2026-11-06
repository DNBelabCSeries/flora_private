# cluster_umis_allbam_v3.py
# Whole-BAM version of Sockeye UMI clustering.
# Processes one contig at a time, writes temporary contig BAMs without indexing,
# then merges -> sorts -> indexes only the final BAM.

import argparse
import collections
import itertools
import logging
import multiprocessing
import os
import pathlib
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pysam
from editdistance import eval as edit_distance
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "bam",
        help="BAM file of alignments with tags for gene (GN), corrected barcode (CB) and uncorrected UMI (UR)",
        type=str,
    )

    parser.add_argument(
        "--output",
        help="Output BAM file with corrected UMI tag (UB) [tagged.sorted.bam]",
        type=str,
        default="tagged.sorted.bam",
    )

    parser.add_argument(
        "-i",
        "--ref_interval",
        help="Size of genomic window (bp) to assign as gene name if no gene assigned by feature annotation [1000]",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--cell_gene_max_reads",
        help="Maximum number of reads to consider for a particular gene + cell barcode combination [20000]",
        type=int,
        default=20000,
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use for per-contig groupby and final sort/index [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    args = parser.parse_args()

    p = pathlib.Path(args.output).resolve()
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.cluster_umis.", dir=str(p.parent))
    args.tempdir_obj = tempdir
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def get_adj_list_directional(umis, counts, threshold=2):
    """
    Identify all UMIs within the Levenshtein distance threshold
    and where the counts of the first UMI is >= (2 * second UMI counts) - 1.
    """
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)

    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


def get_connected_components_adjacency(umis, graph, counts):
    found = set()
    components = []

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)

    return components


def group_directional(clusters, adj_list, counts):
    observed = set()
    groups = []

    for cluster in clusters:
        if len(cluster) == 1:
            groups.append(list(cluster))
            observed.update(cluster)
        else:
            cluster = sorted(cluster, key=lambda x: counts[x], reverse=True)
            temp_cluster = []
            for node in cluster:
                if node not in observed:
                    temp_cluster.append(node)
                    observed.add(node)
            groups.append(temp_cluster)

    return groups


def cluster(counts_dict, threshold=3):
    adj_list = get_adj_list_directional(counts_dict.keys(), counts_dict, threshold)
    clusters = get_connected_components_adjacency(
        counts_dict.keys(), adj_list, counts_dict
    )
    final_umis = [list(x) for x in group_directional(clusters, adj_list, counts_dict)]
    return final_umis


def create_map_to_correct_umi(cluster_list):
    return {y: x[0] for x in cluster_list for y in x}


def correct_umis(umis):
    counts_dict = dict(umis.value_counts())
    umi_map = create_map_to_correct_umi(cluster(counts_dict))
    return umis.replace(umi_map)


def create_region_name(align, args):
    chrom = align.reference_name
    ref_positions = align.get_reference_positions()

    if not ref_positions:
        return "NA"

    start_pos = ref_positions[0]
    end_pos = ref_positions[-1]

    midpoint = int((start_pos + end_pos) / 2)
    interval_start = np.floor(midpoint / args.ref_interval) * args.ref_interval
    interval_end = np.ceil(midpoint / args.ref_interval) * args.ref_interval

    return f"{chrom}_{int(interval_start)}_{int(interval_end)}"


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def launch_pool(func, func_args, procs=1):
    if len(func_args) == 0:
        return []

    if procs <= 1:
        return [func(x) for x in tqdm(func_args, total=len(func_args))]

    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
        raise
    return results


def run_groupby(df):
    df = df.copy()
    df["umi_corr"] = df.groupby(["gene_cell"])["umi_uncorr"].transform(correct_umis)
    return df


def get_bam_info(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    stats = bam.get_index_statistics()
    bam.close()

    n_aligns = int(sum([contig.mapped for contig in stats]))
    chroms = collections.OrderedDict(
        [(contig.contig, contig.mapped) for contig in stats if contig.mapped > 0]
    )
    return n_aligns, chroms


def collect_records_for_contig(input_bam, chrom, args):
    bam = pysam.AlignmentFile(input_bam, "rb")

    records = []
    cell_gene_counter = collections.Counter()

    logger.info(f"[{chrom}] collecting read/gene/barcode/UMI records")
    for align in bam.fetch(contig=chrom):
        if align.is_unmapped:
            continue

        read_id = align.query_name

        for tag in ["UR", "CB", "GN"]:
            assert align.has_tag(tag), f"{tag} tag not found in {read_id}"

        bc_corr = align.get_tag("CB")
        umi_uncorr = align.get_tag("UR")
        gene = align.get_tag("GN")

        if gene == "NA":
            gene = create_region_name(align, args)

        cell_gene_counter[(bc_corr, gene)] += 1
        if args.cell_gene_max_reads == 0 or cell_gene_counter[(bc_corr, gene)] <= args.cell_gene_max_reads:
            records.append((read_id, gene, bc_corr, umi_uncorr))

    bam.close()

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )
    logger.info(f"[{chrom}] collected {df.shape[0]} usable records")
    return df


def compute_corrected_umis_for_contig(df, chrom, args):
    if df.shape[0] == 0:
        return {}, {}

    logger.info(f"[{chrom}] building gene-cell groups")
    df = df.copy()
    df["gene_cell"] = df["gene"] + ":" + df["bc"]
    df = df.set_index("gene_cell")

    gene_cell_unique = list(set(df.index))
    logger.info(f"[{chrom}] total gene-cell groups: {len(gene_cell_unique)}")

    gene_cell_per_chunk = 50
    gene_cell_chunks = list(chunks(gene_cell_unique, gene_cell_per_chunk))

    func_args = []
    for gene_cell_chunk in gene_cell_chunks:
        df_chunk = df.loc[gene_cell_chunk]
        if isinstance(df_chunk, pd.Series):
            df_chunk = df_chunk.to_frame().T
        func_args.append(df_chunk)

    logger.info(f"[{chrom}] clustering UMIs across {len(func_args)} chunks")
    results = launch_pool(run_groupby, func_args, args.threads)

    if len(results) > 0:
        logger.info(f"[{chrom}] concatenating chunk results")
        df = pd.concat(results, axis=0)
    else:
        df = pd.DataFrame(columns=["read_id", "gene", "bc", "umi_uncorr", "umi_corr"])

    logger.info(f"[{chrom}] converting corrected UMIs to dictionaries")
    df = df.drop(["bc", "umi_uncorr"], axis=1).set_index("read_id")
    umis = df["umi_corr"].to_dict()
    genes = df["gene"].to_dict()

    logger.info(f"[{chrom}] corrected reads: {len(umis)}")
    return umis, genes


def write_contig_bam(input_bam, chrom, umis, genes, out_bam_path):
    bam = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(out_bam_path, "wb", template=bam)

    logger.info(f"[{chrom}] writing corrected BAM -> {out_bam_path}")
    for align in bam.fetch(contig=chrom):
        read_id = align.query_name

        if (umis.get(read_id) is not None) and (genes.get(read_id) is not None):
            align.set_tag("UB", umis[read_id], value_type="Z")
            align.set_tag("GN", genes[read_id], value_type="Z")

        bam_out.write(align)

    bam.close()
    bam_out.close()
    logger.info(f"[{chrom}] finished writing temporary BAM")


def merge_contig_bams(contig_bams, output_bam, threads=1):
    if len(contig_bams) == 0:
        raise ValueError("No contig BAMs were produced")

    output_bam = os.path.abspath(output_bam)
    merged_tmp_bam = output_bam.replace(".bam", ".merged.tmp.bam")

    logger.info(f"Merging {len(contig_bams)} contig BAMs into temporary BAM")
    pysam.merge("-f", merged_tmp_bam, *contig_bams)

    logger.info(f"Sorting merged BAM into {output_bam}")
    pysam.sort("-@", str(threads), "-o", output_bam, merged_tmp_bam)

    logger.info(f"Indexing final BAM {output_bam}")
    subprocess.run(["samtools", "index", "-@", str(threads), output_bam], check=True)

    if os.path.exists(merged_tmp_bam):
        os.remove(merged_tmp_bam)


def main(args):
    init_logger(args)

    n_aligns, chroms = get_bam_info(args.bam)
    logger.info(
        f"Input BAM has {n_aligns} mapped alignments across {len(chroms)} contigs"
    )

    contig_bams = []

    for chrom, n_mapped in chroms.items():
        logger.info(f"[{chrom}] start processing ({n_mapped} mapped reads)")

        df = collect_records_for_contig(args.bam, chrom, args)
        umis, genes = compute_corrected_umis_for_contig(df, chrom, args)

        out_bam_path = os.path.join(args.tempdir, f"{chrom}.tagged.bam")
        write_contig_bam(args.bam, chrom, umis, genes, out_bam_path)
        contig_bams.append(out_bam_path)

        logger.info(f"[{chrom}] cleaning up in-memory objects")
        del df
        del umis
        del genes

    logger.info("All contigs finished, starting final merge/sort/index")
    merge_contig_bams(contig_bams, args.output, threads=args.threads)
    logger.info(f"Finished writing {args.output}")


if __name__ == "__main__":
    args = parse_args()
    main(args)
