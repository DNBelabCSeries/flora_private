import importlib.util
import sys
import sqlite3
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ont_features", ROOT / "scripts/ont_features.py")


class OntSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.available = SPEC.origin and Path(SPEC.origin).is_file()
        if cls.available:
            cls.mod = importlib.util.module_from_spec(SPEC)
            sys.modules[SPEC.name] = cls.mod
            SPEC.loader.exec_module(cls.mod)

    def setUp(self):
        self.assertTrue(self.available, "ONT feature backend is not implemented")

    def hit(self, gene="G", transcript="T", score=100, qcov=0.9, tcov=0.9, mapq=60):
        return self.mod.Hit(gene, transcript, score, qcov, tcov, mapq)

    def test_unique_transcript_can_survive_low_genomic_mapq(self):
        result = self.mod.select_features([self.hit()], genome_mapq=1)
        self.assertEqual(result, ("-", "T"))

    def test_unique_transcript_mapq_zero_is_rejected(self):
        self.assertEqual(self.mod.select_features([self.hit(mapq=0)], 60), ("-", "-"))

    def test_gene_threshold_is_inclusive(self):
        self.assertEqual(self.mod.select_features([self.hit()], 30), ("G", "T"))

    def test_score_difference_with_coverage_resolves_transcript(self):
        hits = [self.hit(transcript="T1", qcov=0.4, tcov=0.4), self.hit(transcript="T2", score=99)]
        self.assertEqual(self.mod.select_features(hits, 60), ("G", "T1"))

    def test_score_difference_does_not_bypass_low_coverage(self):
        hits = [self.hit(transcript="T1", qcov=0.39), self.hit(transcript="T2", score=99)]
        self.assertEqual(self.mod.select_features(hits, 60), ("G", "-"))

    def test_identical_best_candidates_remain_transcript_ambiguous(self):
        self.assertEqual(self.mod.select_features([self.hit(transcript="T1"), self.hit(transcript="T2")], 60), ("G", "-"))

    def test_transcript_coverage_tiebreak_must_exceed_point_eight(self):
        for top, expected in [(0.8, "-"), (0.81, "T1")]:
            hits = [self.hit(transcript="T1", tcov=top), self.hit(transcript="T2", tcov=0.7)]
            self.assertEqual(self.mod.select_features(hits, 60), ("G", expected))

    def test_gene_follows_best_ranked_candidate_even_if_transcript_ambiguous(self):
        hits = [self.hit(gene="G1", transcript="T1"), self.hit(gene="G2", transcript="T2")]
        self.assertEqual(self.mod.select_features(hits, 60), ("G1", "-"))

    def test_unique_hit_has_no_multimap_coverage_filter(self):
        self.assertEqual(self.mod.select_features([self.hit(qcov=0.2, tcov=0.2)], 60), ("G", "T"))

    def test_no_reference_candidate_is_unassigned(self):
        self.assertEqual(self.mod.select_features([], 60), ("-", "-"))

    def test_excluded_class_top_hit_is_not_replaced_by_second_candidate(self):
        self.assertEqual(self.mod.select_features([self.hit(transcript="-"), self.hit(score=99)], 60), ("G", "-"))

    def test_reference_join_precedes_candidate_ranking_and_forbidden_class_stays(self):
        with tempfile.TemporaryDirectory() as temp:
            gtf = Path(temp) / "ref.gtf"
            gtf.write_text('chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "G"; transcript_id "T";\n')
            tmap = Path(temp) / "query.tmap"
            tmap.write_text("qry_id\tref_id\tclass_code\nnovel\t-\tu\nknown\tT\t=\nintron\tT\ti\n")
            result = self.mod.load_candidates(tmap, gtf)
            self.assertNotIn("novel", result)
            self.assertEqual(result["intron"], ("G", "-"))
            self.assertEqual(result["known"], ("G", "T"))

    def test_reference_gene_name_can_come_from_gene_or_exon_records(self):
        for feature in ("gene", "exon"):
            for name_first in (True, False):
                with self.subTest(feature=feature, name_first=name_first), tempfile.TemporaryDirectory() as temp:
                    gtf = Path(temp) / "ref.gtf"
                    transcript = 'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T";\n'
                    named = f'chr1\tx\t{feature}\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "GENE_A";\n'
                    gtf.write_text(named + transcript if name_first else transcript + named)
                    self.assertEqual(self.mod.reference_transcripts(gtf), {"T": "GENE_A"})

    def test_reference_uses_gene_id_when_no_gene_name_is_available(self):
        for name in ('', 'gene_name "";', 'gene_name "-";', 'gene_name ".";'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temp:
                gtf = Path(temp) / "ref.gtf"
                gtf.write_text(f'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T"; {name}\n')
                self.assertEqual(self.mod.reference_transcripts(gtf), {"T": "G1"})

    def test_reference_preserves_explicit_transcript_gene_name(self):
        with tempfile.TemporaryDirectory() as temp:
            gtf = Path(temp) / "ref.gtf"
            gtf.write_text('chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "OTHER";\n'
                           'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "GENE_A"; transcript_id "T";\n')
            self.assertEqual(self.mod.reference_transcripts(gtf), {"T": "GENE_A"})

    def test_reference_fallback_prefers_gene_over_exon_and_other_transcript_names(self):
        alias = 'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "ALIAS"; transcript_id "T1";\n'
        gene = 'chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "CANONICAL";\n'
        exon = 'chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "EXON_NAME"; transcript_id "T2";\n'
        unnamed = 'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        for lines in ((alias, exon, gene, unnamed), (gene, unnamed, exon, alias), (exon, unnamed, alias, gene)):
            with self.subTest(lines=lines), tempfile.TemporaryDirectory() as temp:
                gtf = Path(temp) / "ref.gtf"
                gtf.write_text("".join(lines))
                self.assertEqual(self.mod.reference_transcripts(gtf), {"T1": "ALIAS", "T2": "CANONICAL"})

    def test_reference_rejects_transcript_without_gene_identifiers(self):
        with tempfile.TemporaryDirectory() as temp:
            gtf = Path(temp) / "ref.gtf"
            gtf.write_text('chr1\tx\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "T";\n')
            with self.assertRaisesRegex(ValueError, "T"):
                self.mod.reference_transcripts(gtf)

    def test_gene_id_only_candidate_is_written_as_assigned(self):
        import pysam
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            gtf, tmap = root / "ref.gtf", root / "query.tmap"
            gtf.write_text('chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T";\n')
            tmap.write_text("qry_id\tref_id\tclass_code\nassembled\tT\t=\n")
            gene, transcript = self.mod.load_candidates(tmap, gtf)["assembled"]
            selected = self.mod.select_features([self.hit(gene=gene, transcript=transcript)], 60)
            bam = root / "input.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                record = pysam.AlignedSegment(out.header)
                record.query_name, record.query_sequence = "read1", "A" * 50
                record.reference_id, record.reference_start = 0, 0
                record.cigarstring, record.mapping_quality = "50M", 60
                out.write(record)
            with sqlite3.connect(":memory:") as db:
                db.execute("CREATE TABLE features(chrom,read_id,gene,transcript,reason)")
                db.execute("INSERT INTO features VALUES (?,?,?,?,?)", ("chr1", "read1", *selected, "candidate"))
                genes, tx = root / "genes.tsv", root / "tx.tsv"
                self.mod.write_tagged_assignments(db, bam, root / "output.bam", genes, tx, 30, 1)
                self.assertEqual(genes.read_text(), "read1\tAssigned\t60\tG1\n")
                self.assertEqual(tx.read_text(), "read1\tAssigned\t60\tG1\tT\n")

    def test_tagging_uses_primary_chromosome_and_preserves_transcript_only_reads(self):
        import pysam
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            bam = root / "input.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}]}
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for tid, flag in [(0, 0), (1, 2048)]:
                    a = pysam.AlignedSegment(out.header)
                    a.query_name = "read1"
                    a.query_sequence = "A" * 50
                    a.flag, a.reference_id, a.reference_start = flag, tid, 0
                    a.cigarstring, a.mapping_quality = "50M", 1
                    a.set_tags([("CB", "cell1"), ("UR", "AAAA")])
                    out.write(a)
            db = sqlite3.connect(":memory:")
            db.execute("CREATE TABLE features(chrom,read_id,gene,transcript,reason)")
            db.executemany("INSERT INTO features VALUES (?,?,?,?,?)", [("chr1", "read1", "-", "T1", "candidate"), ("chr2", "read1", "G2", "T2", "candidate")])
            output = root / "output.bam"
            genes, tx = root / "genes.tsv", root / "tx.tsv"
            self.mod.write_tagged_assignments(db, bam, output, genes, tx, 30, 1)
            with pysam.AlignmentFile(output) as result:
                records = list(result)
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].get_tag("GN"), "NA")
            self.assertEqual(records[0].get_tag("TX"), "T1")
            self.assertIn("Unassigned_mapq", genes.read_text())
            self.assertIn("\tAssigned\t1\tNA\tT1", tx.read_text())
            db.close()


if __name__ == "__main__":
    unittest.main()
