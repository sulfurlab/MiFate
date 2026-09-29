"""Frozen MiFate deployment: score KO matrices for SEC and DIS."""
from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAIN_COLUMNS = [
    "Genome_ID", "SEC_score", "SEC_prediction", "SEC_domain",
    "DIS_score", "DIS_prediction", "DIS_domain",
]


def read_table(path):
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path)
    return pd.read_csv(path, sep=None, engine="python")


def ids(df, name):
    if len(df.columns) < 2:
        raise ValueError(f"{name}: expected ID column followed by data columns")
    df = df.copy()
    first = df.columns[0]
    if df[first].isna().any():
        raise ValueError(f"{name}: missing MAG ID")
    df[first] = df[first].astype(str).str.strip()
    if df[first].eq("").any() or df[first].duplicated().any():
        raise ValueError(f"{name}: empty or duplicate MAG ID")
    return df.set_index(first).rename_axis("Genome_ID")


def quality_gate(q, index):
    lookup = {str(c).strip().lower(): c for c in q.columns}
    if not {"completeness", "contamination"} <= lookup.keys():
        raise ValueError("Quality file requires Completeness and Contamination columns (percent units 0-100).")
    out = pd.DataFrame(index=index)
    for column in ("Completeness", "Contamination"):
        out[column] = pd.to_numeric(q[lookup[column.lower()]], errors="coerce").reindex(index)
    finite = np.isfinite(out.Completeness) & np.isfinite(out.Contamination)
    valid = finite & out.Completeness.between(0, 100) & out.Contamination.ge(0)
    out["QC_pass"] = valid & (out.Completeness >= 50) & (out.Contamination < 10)
    out["QC_reason"] = np.where(
        ~valid, "Missing or invalid quality values",
        np.where(out.QC_pass, "", "Requires completeness ≥50% and contamination <10%"),
    )
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ko", type=Path, required=True)
    qc = p.add_mutually_exclusive_group(required=True)
    qc.add_argument("--quality", type=Path)
    qc.add_argument("--assume-qc-pass", action="store_true",
                    help="Use only after confirming every input genome meets completeness ≥50%% and contamination <10%%")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--artifacts", type=Path, default=ROOT / "models/applicability_domain")
    p.add_argument("--rules", type=Path, default=ROOT / "models/deployment_rules.json")
    a = p.parse_args()

    required = [a.rules]
    for stage in ("SEC", "DIS"):
        required.extend(a.artifacts / stage / name for name in ("model.json", "reference.json", "reference.npz"))
    missing_artifacts = [str(path) for path in required if not path.is_file()]
    if missing_artifacts:
        p.error("Missing frozen MiFate artifacts: " + ", ".join(missing_artifacts))
    try:
        import applicability_domain as ad
    except ModuleNotFoundError as exc:
        p.error(f"Missing prediction dependency: {exc}. Install requirements.txt first.")

    k = ids(read_table(a.ko), "KO matrix")
    if a.assume_qc_pass:
        result = pd.DataFrame(index=k.index)
        result["Completeness"] = np.nan
        result["Contamination"] = np.nan
        result["QC_pass"] = True
        result["QC_reason"] = ""
        result["QC_source"] = "user assertion: all input genomes pass QC"
    else:
        q = ids(read_table(a.quality), "Quality")
        result = quality_gate(q, k.index)
        result["QC_source"] = "provided quality table"
    rules = json.loads(a.rules.read_text(encoding="utf-8"))
    audit = {}
    a.output.mkdir(parents=True, exist_ok=True)

    for stage in ("SEC", "DIS"):
        bundle = a.artifacts / stage
        _, features, _, ref, dist = ad.load(bundle)
        r = rules[stage]
        # Existing frozen rule files call this field persistence_min.
        retention_min = r.get("retention_min", r.get("persistence_min"))
        if retention_min is None:
            raise ValueError(f"{stage}: missing retention_min/persistence_min rule")

        f = k.copy()
        f.columns = [str(c).strip().removeprefix("KEGG_") for c in f.columns]
        names = [str(c).removeprefix("KEGG_") for c in features]
        if not f.columns.is_unique:
            raise ValueError("Duplicate normalized KO column names")
        missing = sorted(set(names) - set(f.columns))
        extra = sorted(set(f.columns) - set(names))
        x = f.reindex(columns=names, fill_value=0).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        good = np.isfinite(x).all(axis=1) & np.isin(x, [0, 1]).all(axis=1) & (np.nansum(x, axis=1) > 0)
        eligible = result.QC_pass.to_numpy() & good

        out = pd.DataFrame(index=k.index)
        out["score"] = np.nan
        out["domain"] = pd.Series(pd.NA, index=k.index, dtype="object")
        out["domain_percentile"] = np.nan
        out["domain_distance"] = np.nan
        out["prediction"] = "Not predicted"
        out["reason"] = result.QC_reason
        out.loc[result.QC_pass & ~good, "reason"] = "Invalid binary KO values or no present model KOs"

        if eligible.any():
            xx = x[eligible].astype(np.uint8)
            dd = ad.measure(xx, ref)[:, 1]
            pct = 100 * np.searchsorted(np.sort(dist[:, 1]), dd, side="right") / len(dist)
            score = ad.model_probability(bundle / "model.json", xx, features)
            outside = pct > r["max_domain_percentile"]
            labels = np.where(
                outside, "Uncertain",
                np.where(score <= r["attenuation_max"], "Relative attenuation",
                         np.where(score >= retention_min, "Relative retention", "Uncertain")),
            )
            out.loc[eligible, "score"] = score
            out.loc[eligible, "domain"] = np.where(outside, "Out-of-domain", "In-domain")
            out.loc[eligible, "domain_percentile"] = pct
            out.loc[eligible, "domain_distance"] = dd
            out.loc[eligible, "prediction"] = labels
            out.loc[eligible, "reason"] = np.where(
                outside, "Outside applicability domain",
                np.where(labels == "Uncertain", "Score between class thresholds", ""),
            )

        result = result.join(out.add_prefix(stage + "_"))
        audit[stage] = dict(
            model_sha256=ad.sha(bundle / "model.json"),
            reference_sha256=ad.sha(bundle / "reference.npz"),
            rules=r,
            missing_KOs_zero_filled=missing,
            extra_KOs_ignored=extra,
            counts=out.prediction.value_counts().to_dict(),
        )

    detailed = result.reset_index()
    predicted = detailed.SEC_prediction.ne("Not predicted") & detailed.DIS_prediction.ne("Not predicted")
    detailed.loc[predicted, MAIN_COLUMNS].to_csv(
        a.output / "MiFate_predictions.csv", index=False, encoding="utf-8-sig"
    )
    detailed.to_csv(a.output / "MiFate_predictions_detailed.csv", index=False, encoding="utf-8-sig")
    detailed.loc[~predicted, [
        "Genome_ID", "Completeness", "Contamination", "QC_pass", "QC_source", "QC_reason",
        "SEC_reason", "DIS_reason",
    ]].to_csv(a.output / "excluded_genomes.csv", index=False, encoding="utf-8-sig")

    (a.output / "run_audit.json").write_text(json.dumps(dict(
        ko_file=str(a.ko.resolve()), quality_file=str(a.quality.resolve()) if a.quality else None,
        ko_sha256=ad.sha(a.ko), quality_sha256=ad.sha(a.quality) if a.quality else None,
        QC_source="user assertion: all input genomes pass QC" if a.assume_qc_pass else "provided quality table",
        N=len(k), QC_pass=int(result.QC_pass.sum()),
        n_predicted=int(predicted.sum()), n_excluded=int((~predicted).sum()),
        stages=audit,
    ), ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Completed: {len(k)} genomes; predicted: {predicted.sum()}; excluded: {(~predicted).sum()}. Output: {a.output.resolve()}")
    for stage in ("SEC", "DIS"):
        print(stage, audit[stage]["counts"])


if __name__ == "__main__":
    main()
