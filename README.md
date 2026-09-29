# MiFate: genome-informed treatment-response predictions

MiFate takes a microbial genome FASTA file (or a batch), performs protein and
KO annotation, and returns separate predictions for secondary effluent (SEC)
and post-disinfection effluent (DIS). The categories are **Relative
attenuation**, **Relative retention**, and **Uncertain**. The output describes
general response tendencies learned from community-normalized observations.

**Model version:** this repository includes the **2026-09-29 retrained model**.
It is not the earlier frozen manuscript model. See [MODEL_CARD.md](MODEL_CARD.md)
for metrics and interpretation. No original MAG sequences, WHO input genomes,
or large KOfam databases are distributed here.

## Repository contents

```text
bin/run_mifate_full.sh                   optional Bash wrapper
src/genome_deployment_pipeline.py        genome FASTA → KO matrix → prediction
src/mifate_predict.py                    existing binary KO matrix → prediction
src/applicability_domain.py              model load, Jaccard distance, XGBoost
models/deployment_rules.json             version-specific score cutoffs
models/applicability_domain/{SEC,DIS}/   model.json, reference.npz, reference.json
examples/genome_quality.csv              input template
tests/test_deployment.py                 KO threshold and model integration checks
validation/                              training and confidence-interval summaries
```

## Requirements

- Python 3.10 or newer: `python -m pip install -r requirements.txt`
- [Prodigal](https://github.com/hyattpd/Prodigal), available as `prodigal`
- [KOfamScan](https://github.com/takaram/kofam_scan), available as
  `exec_annotation`; HMMER and its other dependencies must also be installed
- KOfam profile HMMs (`profiles/` or a `.hal` file) **and** the corresponding
  `ko_list` file containing the KO-specific score thresholds

The profile HMM path and `ko_list` path are different inputs. The KOfamScan
installation and database are supplied by the operator; they are not bundled.
Confirm that their release and the gene-prediction procedure are compatible
with the training KO matrix before relying on new-genome comparisons. The
training matrix contains binary KO calls but no raw KOfamScan scores, so this
compatibility cannot yet be certified from the supplied files alone.

## Genome FASTA input

Use one `.fna`, `.fa`, or `.fasta` per genome. A directory or ZIP archive of
such files is also accepted. Filename stems must be unique. Supply a CSV/TSV
quality table with matching IDs and percent units:

```csv
genome_id,completeness,contamination
MAG_001,82.4,2.1
MAG_002,50.0,9.9
```

The gate is **completeness ≥50% and contamination <10%**. Genomes failing it
are reported as excluded, without running annotation or classification.

Run from the repository root:

```bash
python src/genome_deployment_pipeline.py \
  --genomes /data/new_genomes \
  --quality /data/genome_quality.csv \
  --kofam-profile /opt/kofam/profiles \
  --kofam-ko-list /opt/kofam/ko_list \
  --output /data/mifate_run_001 \
  --threads 8 \
  --keep-intermediates
```

Or use the Bash wrapper on Linux/macOS:

```bash
bash bin/run_mifate_full.sh \
  -i /data/new_genomes -q /data/genome_quality.csv \
  -d /opt/kofam/profiles -k /opt/kofam/ko_list \
  -o /data/mifate_run_001 -t 8
```

The output directory must be new or empty. `--keep-intermediates` retains
Prodigal proteins, KOfamScan `detail-tsv` output, and the aligned KO matrix
for inspection. Without it, the work directory is removed after a successful
run. `pipeline.log` and `run_manifest.json` retain annotation settings and
accepted-KO counts. This package requests `detail-tsv` and **uses only `*`
marked KOfamScan hits**. Unmarked hits are not considered present KOs. It
does not change KOfamScan's default E-value or threshold-scale options.

## Precomputed KO matrix input

If each input genome has already passed QC and has a **binary** KO profile:

```bash
python src/mifate_predict.py \
  --ko /data/new_KO_matrix.csv \
  --quality /data/genome_quality.csv \
  --output /data/mifate_matrix_run_001
```

Use `--assume-qc-pass` instead of `--quality` only when the QC status of every
input genome is independently established. The first CSV column is the genome
ID; all remaining columns are KO IDs with 0/1 values. Missing model KOs are
zero-filled, extra KOs discarded, and columns aligned to the saved model
feature order. **This entry point cannot infer whether a supplied matrix was
built using starred KOfamScan hits**; retain the raw annotation provenance.

## Outputs

```text
mifate_run_001/
  predictions/
    MiFate_predictions.csv             7-column main output, predicted genomes
    MiFate_predictions_detailed.csv    QC, domain distances and reasons
    excluded_genomes.csv               QC/KO failures
    run_audit.json                     model hashes, rule values, KO alignment
  qc_summary.csv
  run_manifest.json
  pipeline.log
  work/                                  only with --keep-intermediates
```

The main columns are `Genome_ID, SEC_score, SEC_prediction, SEC_domain,
DIS_score, DIS_prediction, DIS_domain`. For each stage, genomes exceeding the
training 99th-percentile mean five-nearest-neighbor Jaccard distance are
`Uncertain`. In-domain scores between that stage's two fixed class cutoffs are
also `Uncertain`. There is no retraining during prediction.

## Checks and GitHub publication

```bash
python -m unittest discover -s tests -v
bash -n bin/run_mifate_full.sh
sha256sum -c MANIFEST_SHA256.txt
```

The tests use a small mock KOfamScan output to check the `*` rule and the real
bundled model to check inference. They do not replace a real-data annotation
check with your installed Prodigal/KOfamScan and database.

This directory can be copied into a new GitHub repository and pushed after
reviewing the model artifact release rights. The included reference files
contain training feature profiles and training genome identifiers. Add an
appropriate project license and manuscript citation before making the
repository public. Do not commit input genomes or generated `work/` results.

After creating an empty GitHub repository, push this folder with:

```bash
git init
git add .
git commit -m "Add MiFate retrained deployment pipeline"
git branch -M main
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
```

Replace `OWNER/REPOSITORY` with your repository address. GitHub Actions runs
the tests and checksum check on pushes and pull requests. Predictions from
large genome batches still run on a workstation or compute cluster supplied
with Prodigal, KOfamScan, and the KOfam database.
