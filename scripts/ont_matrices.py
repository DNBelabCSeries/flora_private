"""Export Flora counts as MEX, with wf-single-cell v3.3.0 preprocessing."""

import argparse
import csv
import gzip
import json
from pathlib import Path


def read_counts(path):
    import numpy as np
    from scipy import sparse

    features, row_idx, col_idx, values = [], [], [], []
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        cells = header[1:]
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"Malformed matrix row in {path}: {row[0]}")
            features.append(row[0])
            for col, value in enumerate(row[1:]):
                count = int(value)
                if count < 0:
                    raise ValueError("Expression counts must be nonnegative")
                if count:
                    row_idx.append(len(features) - 1)
                    col_idx.append(col)
                    values.append(count)
    matrix = sparse.csr_matrix((np.array(values, dtype=np.int64), (row_idx, col_idx)), shape=(len(features), len(cells)))
    return matrix, features, cells


def process_matrix(matrix, features, cells, min_features=200, min_cells=3, max_mito=20, norm_count=10000,
                   mito_prefixes=("MT-", "mt-")):
    import numpy as np
    from scipy import sparse

    cell_mask = np.asarray(matrix.getnnz(axis=0)) >= min_features
    feature_mask = np.asarray(matrix.getnnz(axis=1)) >= min_cells
    selected_features = [name for name, keep in zip(features, feature_mask) if keep]
    selected_cells = [name for name, keep in zip(cells, cell_mask) if keep]
    filtered = matrix[feature_mask, :][:, cell_mask].astype(float).tocsr()
    if filtered.shape[1]:
        totals = np.asarray(filtered.sum(axis=0)).ravel()
        mito = [any(prefix in name for prefix in mito_prefixes) for name in selected_features]
        mito_counts = np.asarray(filtered[mito, :].sum(axis=0)).ravel()
        keep = (totals > 0) & (mito_counts <= (max_mito / 100) * totals)
        selected_cells = [name for name, flag in zip(selected_cells, keep) if flag]
        filtered = filtered[:, keep]
        totals = totals[keep]
        if len(totals):
            filtered = (filtered @ sparse.diags(norm_count / totals)).tocsr()
            filtered.data = np.log1p(filtered.data) / np.log(10)
    return filtered, selected_features, selected_cells


def write_mex(directory, matrix, features, cells, feature_type):
    from scipy.io import mmwrite

    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, matrix, symmetry="general")
    with gzip.open(directory / "features.tsv.gz", "wt") as handle:
        for feature in features:
            handle.write(f"{feature}\t{feature}\t{feature_type}\n")
    with gzip.open(directory / "barcodes.tsv.gz", "wt") as handle:
        handle.writelines(f"{cell}\n" for cell in cells)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-tsv", required=True, type=Path)
    parser.add_argument("--transcript-tsv", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    stats = {"reference": "wf-single-cell v3.3.0", "reference_commit": "d46f92121c79c28a02105cf2444029e91dc4b609",
             "reference_withdrawn": True, "normalization": "log10(1 + counts * 10000 / cell_total)",
             "min_features": 200, "min_cells": 3, "max_mito_percent": 20,
             "mito_identifiers": ["MT-", "mt-"], "report_population": "Flora original cell calls; processed MEX has additional QC"}
    for feature, path in [("gene", args.gene_tsv), ("transcript", args.transcript_tsv)]:
        if path is None:
            continue
        matrix, features, cells = read_counts(path)
        kind = "Gene Expression" if feature == "gene" else "Transcript Expression"
        write_mex(args.output_dir / f"{args.sample}.{feature}_raw_feature_bc_matrix", matrix, features, cells, kind)
        processed, pfeatures, pcells = process_matrix(matrix, features, cells)
        write_mex(args.output_dir / f"{args.sample}.{feature}_processed_feature_bc_matrix", processed, pfeatures, pcells, kind)
        stats[feature] = {"raw_features": len(features), "raw_cells": len(cells),
                          "processed_features": len(pfeatures), "processed_cells": len(pcells)}
    (args.output_dir / f"{args.sample}.matrix_processing.json").write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    main()
