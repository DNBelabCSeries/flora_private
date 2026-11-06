"""Small real-tool smoke test; runs when the assembly toolchain is on PATH."""

import argparse
import importlib.util
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest


TOOLS = ("samtools", "minimap2", "stringtie", "gffread", "gffcompare")
ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(all(shutil.which(tool) for tool in TOOLS), "requires ONT assembly toolchain")
class ExternalSmokeTests(unittest.TestCase):
    def test_assembly_remapping_and_primary_counting(self):
        import pysam

        spec = importlib.util.spec_from_file_location("ont_features_external", ROOT / "scripts/ont_features.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            genome = root / "genome.fa"
            rng = random.Random(42)
            dna = "".join(rng.choices("ACGT", k=5000))
            genome.write_text(f">chr1\n{dna}\n>chr2\n{dna[::-1]}\n")
            gtf = root / "genes.gtf"
            lines = []
            for chrom, gene, tx in [("chr1", "G1", "T1"), ("chr2", "G2", "T2")]:
                attr = f'gene_id "{gene}"; gene_name "{gene}"; transcript_id "{tx}";'
                for feature, start, end in [("transcript", 1001, 2000), ("exon", 1001, 1300), ("exon", 1601, 2000)]:
                    lines.append(f"{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\t{attr}\n")
            gtf.write_text("".join(lines))
            bam = root / "reads.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 5000}, {"SN": "chr2", "LN": 5000}]}
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for i in range(20):
                    rec = pysam.AlignedSegment(out.header)
                    rec.query_name = f"read{i}"
                    rec.query_sequence = dna[1000:1300] + dna[1600:2000]
                    rec.query_qualities = pysam.qualitystring_to_array("I" * 700)
                    rec.reference_id, rec.reference_start, rec.flag = 0, 1000, 0
                    rec.mapping_quality = 1 if i == 0 else 60
                    rec.cigarstring = "300M300N400M"
                    rec.set_tags([("CB", f"cell{i % 4}"), ("UR", "ACGTACGTAA"), ("XS", "+", "A")])
                    out.write(rec)
                # An isolated alignment outside chr2 annotation must not become a gene.
                rec = pysam.AlignedSegment(out.header)
                rec.query_name = "no_feature"
                rec.query_sequence = dna[::-1][3000:3700]
                rec.query_qualities = pysam.qualitystring_to_array("I" * 700)
                rec.reference_id, rec.reference_start, rec.flag = 1, 3000, 0
                rec.mapping_quality, rec.cigarstring = 60, "700M"
                rec.set_tags([("CB", "cell0"), ("UR", "CCCCCCCCCC")])
                out.write(rec)
            pysam.index(str(bam))
            args = argparse.Namespace(bam=bam, tagged_bam=bam, genome=genome, gtf=gtf,
                                      output_dir=root / "assembly", output_bam=root / "tagged.bam",
                                      gene_assigns=root / "gene.tsv", transcript_assigns=root / "tx.tsv",
                                      sample="test", threads=4, mapq=30, stringtie_opts="-c 2", save_intermediate=False)
            mod.run(args)
            genes = [line.split("\t") for line in args.gene_assigns.read_text().splitlines()]
            txs = [line.split("\t") for line in args.transcript_assigns.read_text().splitlines()]
            self.assertEqual(len(genes), 21)
            self.assertEqual(sum(row[1] == "Assigned" for row in genes), 19)
            self.assertEqual(sum(row[1] == "Assigned" for row in txs), 20)
            self.assertFalse(Path(str(genome) + ".fai").exists())
            self.assertFalse(list(args.output_dir.glob("ont-work-*")))
            self.assertTrue((args.output_dir / "test.transcriptome.fa.gz").is_file())
            # A stale link from a previous run must never select an old reference.
            (args.output_dir / "reference.fa").symlink_to(root / "missing-old-genome.fa")
            args.save_intermediate = True
            mod.run(args)
            serial_genes = args.gene_assigns.read_bytes()
            serial_txs = args.transcript_assigns.read_bytes()
            with pysam.AlignmentFile(args.output_bam) as handle:
                serial_records = [r.to_string() for r in handle]
            args.jobs = 2
            args.threads = 16
            args.temp_dir = root / "local_scratch"
            mod.run(args)
            self.assertEqual(args.gene_assigns.read_bytes(), serial_genes)
            self.assertEqual(args.transcript_assigns.read_bytes(), serial_txs)
            with pysam.AlignmentFile(args.output_bam) as handle:
                self.assertEqual([r.to_string() for r in handle], serial_records)
            self.assertFalse(list(args.temp_dir.iterdir()))
            saved = list((args.output_dir / "intermediates").glob("ont-work-*/chrom_00000/reference.gtf"))
            self.assertEqual(len(saved), 2)
            stage_dir = os.environ.get("FLORA_TEST_BIN_DIR")
            if stage_dir:
                corrected = root / "corrected.bam"
                table, gene_matrix, tx_matrix = root / "table.tsv", root / "gene_matrix.tsv", root / "tx_matrix.tsv"
                subprocess.run([str(Path(stage_dir) / "cluster_umis_allbam"), "--cell_gene_max_reads", "0",
                                "--output", str(corrected), str(args.output_bam)], check=True)
                subprocess.run([str(Path(stage_dir) / "cell_umi_gene_table"), "--output", str(table),
                                "--gene-expression-output", str(gene_matrix), str(corrected)], check=True)
                subprocess.run([str(Path(stage_dir) / "isoform_expression"), "--output", str(tx_matrix),
                                str(corrected), str(args.transcript_assigns)], check=True)
                self.assertEqual(len(table.read_text().splitlines()), 22)
                self.assertEqual(gene_matrix.read_text().splitlines()[1], "G1\t1\t1\t1\t1")
                self.assertEqual(tx_matrix.read_text().splitlines()[1], "T1\t1\t1\t1\t1")


if __name__ == "__main__":
    unittest.main()
