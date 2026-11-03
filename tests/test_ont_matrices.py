import importlib.util
import gzip
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("ont_matrices", ROOT / "scripts/ont_matrices.py")
        cls.mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.mod)

    def test_filters_are_simultaneous_then_normalized_log10(self):
        script = ROOT / "scripts/ont_matrices.py"
        self.assertTrue(script.exists(), "ONT sparse matrix exporter is missing")
        spec = importlib.util.spec_from_file_location("ont_matrices", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        import numpy as np
        from scipy import sparse
        raw = sparse.csr_matrix([[1, 1, 0], [1, 0, 0], [0, 1, 5]])
        matrix, features, cells = mod.process_matrix(raw, ["G1", "G2", "MT-G"], ["C1", "C2", "C3"],
                                                     min_features=2, min_cells=2, max_mito=40, norm_count=100)
        # Initial masks retain G1,MT-G and C1,C2. C2 then fails mito filtering.
        self.assertEqual(features, ["G1", "MT-G"])
        self.assertEqual(cells, ["C1"])
        np.testing.assert_allclose(matrix.toarray(), [[np.log10(101)], [0]])

    def test_no_retained_features_removes_zero_total_cells_and_barcodes(self):
        from scipy import sparse
        from scipy.io import mmread
        raw = sparse.csr_matrix([[1]] * 200)
        matrix, features, cells = self.mod.process_matrix(raw, [f"G{i}" for i in range(200)], ["C1"])
        self.assertEqual(features, [])
        self.assertEqual(cells, [])
        self.assertEqual(matrix.shape, (0, 0))
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp) / "processed"
            self.mod.write_mex(directory, matrix, features, cells, "Gene Expression")
            with gzip.open(directory / "barcodes.tsv.gz", "rt") as handle:
                self.assertEqual(handle.read(), "")
            with gzip.open(directory / "matrix.mtx.gz", "rb") as handle:
                self.assertEqual(mmread(handle).shape, (0, 0))

    def test_zero_total_cell_removal_is_independent_of_retained_feature_count(self):
        from scipy import sparse
        for extra_feature in (False, True):
            with self.subTest(extra_feature=extra_feature):
                rows = [[1, 0, 0, 0] for _ in range(200)]
                if extra_feature:
                    rows.append([0, 1, 1, 1])
                matrix, features, cells = self.mod.process_matrix(
                    sparse.csr_matrix(rows), [f"G{i}" for i in range(len(rows))], ["C1", "C2", "C3", "C4"])
                self.assertEqual(cells, [])
                self.assertEqual(matrix.shape, (int(extra_feature), 0))

    def test_initially_empty_matrices_remain_empty(self):
        from scipy import sparse
        for shape in ((0, 0), (0, 1), (1, 0)):
            with self.subTest(shape=shape):
                matrix, features, cells = self.mod.process_matrix(
                    sparse.csr_matrix(shape), ["G1"][:shape[0]], ["C1"][:shape[1]],
                    min_features=0, min_cells=0)
                self.assertEqual(cells, [])
                self.assertEqual(matrix.shape, (shape[0], 0))


if __name__ == "__main__":
    unittest.main()
