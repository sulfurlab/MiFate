"""Frozen MiFate model loading, five-neighbor Jaccard distance, and scoring.

This module is specific to the retrained sparse-CSR XGBoost bundles.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy import sparse
import xgboost as xgb


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def measure(query, reference, exclude=None):
    """Return 1-NN Jaccard, mean 5-NN Jaccard, and 1-NN cosine distances."""
    if len(reference) < 6:
        raise ValueError("At least six training genomes are required.")
    ref = sparse.csr_matrix(reference, dtype=np.float32)
    rs = np.asarray(ref.sum(axis=1)).ravel()
    out = []
    for start in range(0, len(query), 128):
        q = sparse.csr_matrix(query[start:start + 128], dtype=np.float32)
        qs = np.asarray(q.sum(axis=1)).ravel()
        intersection = (q @ ref.T).toarray()
        union = qs[:, None] + rs[None, :] - intersection
        jaccard = 1 - np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0)
        denominator = np.sqrt(qs[:, None] * rs[None, :])
        cosine = 1 - np.divide(intersection, denominator,
                               out=np.zeros_like(intersection), where=denominator > 0)
        jaccard[union == 0] = 0
        cosine[(qs[:, None] == 0) & (rs[None, :] == 0)] = 0
        if exclude is not None:
            for i, row in enumerate(exclude[start:start + 128]):
                if row >= 0:
                    jaccard[i, row] = np.inf
                    cosine[i, row] = np.inf
        five = np.partition(jaccard, 4, axis=1)[:, :5]
        out.append(np.column_stack((five.min(axis=1), five.mean(axis=1), cosine.min(axis=1))))
    return np.vstack(out) if out else np.empty((0, 3))


def model_probability(modelpath, x, features):
    """Score using the same sparse zero-as-absent representation as training."""
    model = xgb.Booster()
    model.load_model(modelpath)
    matrix = sparse.csr_matrix(np.asarray(x, dtype=np.float32))
    return model.predict(xgb.DMatrix(matrix, feature_names=list(features)))


def load(bundle):
    bundle = Path(bundle)
    cfg = json.loads((bundle / "reference.json").read_text(encoding="utf-8"))
    with np.load(bundle / "reference.npz", allow_pickle=False) as z:
        features = z["features"].tolist()
        training_ids = z["training_ids"].tolist()
        x = z["X"].copy()
        distances = z["training_distances"].copy()
    if sha(bundle / "model.json") != cfg["model_sha256"]:
        raise ValueError("Frozen model checksum mismatch.")
    if sha(bundle / "reference.npz") != cfg["reference_sha256"]:
        raise ValueError("Domain reference checksum mismatch.")
    if x.shape[1] != len(features) or cfg["n_features"] != len(features):
        raise ValueError("Feature order/width mismatch in frozen bundle.")
    return cfg, features, training_ids, x, distances
