"""Version provenance and optional comparison against the original ONT source.

Set ONT_V330_REFERENCE_DIR to an unpacked d46f921 source tree to run the
original assignment script as an oracle. ONT source is not bundled in Flora.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import types
import unittest


ROOT = Path(__file__).resolve().parents[1]
COMMIT = "d46f92121c79c28a02105cf2444029e91dc4b609"


def load_flora():
    spec = importlib.util.spec_from_file_location("flora_ont_v330", ROOT / "scripts/ont_features.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ProvenanceTests(unittest.TestCase):
    def test_backend_has_pinned_withdrawn_reference(self):
        module = load_flora()
        self.assertEqual(module.REFERENCE_VERSION, "wf-single-cell v3.3.0")
        self.assertEqual(module.REFERENCE_COMMIT, COMMIT)
        self.assertTrue(module.REFERENCE_WITHDRAWN)

    def test_build_and_runtime_environments_pin_stringtie_222(self):
        for name in ("environment.yml", "environment.runtime.yml"):
            self.assertIn("- stringtie=2.2.2\n", (ROOT / name).read_text())

    def test_both_entrypoints_identify_330(self):
        for name in ("run_all.sh", "run_all_mixed_species.sh"):
            text = (ROOT / name).read_text()
            self.assertIn("wf-single-cell v3.3.0", text)
            self.assertNotIn("wf-single-cell v3.3.4", text)
            self.assertIn(COMMIT, text)

    def test_legacy_parameter_table_does_not_claim_ont_provenance(self):
        for name in ("run_all.sh", "run_all_mixed_species.sh"):
            text = (ROOT / name).read_text()
            function = "write_parameters_tsv() {" + text.split("write_parameters_tsv() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"
            with tempfile.TemporaryDirectory() as temp:
                for backend in ("ont", "legacy"):
                    output = Path(temp) / f"{backend}.tsv"
                    subprocess.run(["bash", "-c", function + '\nFEATURE_BACKEND="$1"\nwrite_parameters_tsv "$2"\n',
                                    "bash", backend, str(output)], check=True, env={"PATH": "/usr/bin:/bin"})
                    with output.open() as handle:
                        rows = dict(list(csv.reader(handle, delimiter="\t"))[1:])
                    if backend == "ont":
                        self.assertEqual(rows["feature_reference_commit"], COMMIT)
                    else:
                        self.assertEqual(rows["feature_reference"], "Flora GTF-overlap")
                        self.assertNotIn("feature_reference_commit", rows)
                        self.assertNotIn("feature_reference_withdrawn", rows)


@unittest.skipUnless(os.environ.get("ONT_V330_REFERENCE_DIR"), "requires original ONT 3.3.0 source and oracle dependencies")
class OriginalSourceParityTests(unittest.TestCase):
    def test_500_read_assignments_match_original_source(self):
        import pysam

        source = Path(os.environ["ONT_V330_REFERENCE_DIR"])
        manifest = json.loads((ROOT / "docs/ont-v330-reference.json").read_text())
        for relative, expected in manifest["sha256"].items():
            self.assertEqual(hashlib.sha256((source / relative).read_bytes()).hexdigest(), expected)
        # Only logger/parser helpers are replaced; the actual ONT parsing,
        # joining, ranking and selection functions execute unchanged.
        package = types.ModuleType("_ont330_oracle")
        package.__path__ = [str(source / "bin/workflow_glue")]
        util = types.ModuleType("_ont330_oracle.util")
        util.get_named_logger = logging.getLogger
        util.wf_parser = argparse.ArgumentParser
        sys.modules[package.__name__] = package
        sys.modules[util.__name__] = util
        spec = importlib.util.spec_from_file_location("_ont330_oracle.assign_features", source / "bin/workflow_glue/assign_features.py")
        oracle = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = oracle
        spec.loader.exec_module(oracle)
        flora = load_flora()
        rng = random.Random(330)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            bam = root / "transcript.bam"
            gtf, tmap, tags, output = [root / name for name in ("ref.gtf", "query.tmap", "tags.tsv", "assigned.tsv")]
            gtf.write_text('chr1\ttest\ttranscript\t1\t2000\t.\t+\t.\tgene_id "g1"; gene_name "G1"; transcript_id "T1";\n'
                           'chr1\ttest\ttranscript\t1\t2000\t.\t+\t.\tgene_id "g2"; gene_name "G2"; transcript_id "T2";\n')
            tmap.write_text("qry_id\tref_id\tclass_code\nm1\tT1\t=\nm2\tT2\t=\nmi\tT1\ti\nnovel\t-\tu\n")
            lengths = {"m1": 1000, "m2": 800, "mi": 900, "novel": 1200}
            names = list(lengths)
            header = {"HD": {"SO": "queryname"}, "SQ": [{"SN": n, "LN": size} for n, size in lengths.items()]}
            candidates = flora.load_candidates(tmap, gtf)
            expected = {}
            with pysam.AlignmentFile(bam, "wb", header=header) as dest, tags.open("w") as tag_out:
                tag_out.write("read_id\tmapq\n")
                for i in range(500):
                    name = f"read{i:05d}"
                    genome_mapq = rng.choice([0, 1, 29, 30, 60])
                    tag_out.write(f"{name}\t{genome_mapq}\n")
                    hits = []
                    for j, target in enumerate(rng.sample(names, rng.randint(1, 4))):
                        rec = pysam.AlignedSegment(dest.header)
                        rec.query_name, rec.query_sequence = name, "A" * 1000
                        rec.flag = 0 if j == 0 else 256
                        rec.reference_id, rec.reference_start = names.index(target), 0
                        aligned = rng.randint(100, min(950, lengths[target]))
                        rec.cigarstring = f"{aligned}M{1000 - aligned}S"
                        rec.mapping_quality = rng.choice([0, 1, 30, 60])
                        rec.set_tag("AS", rng.randrange(1000000))
                        dest.write(rec)
                        hit = flora.make_hit(rec, candidates, lengths)
                        if hit is not None:
                            hits.append(hit)
                    expected[name] = flora.select_features(hits, genome_mapq)
            args = argparse.Namespace(transcriptome_bam=bam, gffcompare_tmap=tmap, gtf=gtf,
                                      tags=tags, output=output, chunksize=137,
                                      min_tr_coverage=0.4, min_read_coverage=0.4, min_mapq=30)
            oracle.main(args)
            with output.open() as handle:
                reader = csv.reader(handle, delimiter="\t")
                next(reader)
                observed = {row[0]: tuple(row[1:]) for row in reader}
            for name, labels in expected.items():
                self.assertEqual(observed.get(name, ("-", "-")), labels, name)


if __name__ == "__main__":
    unittest.main()
