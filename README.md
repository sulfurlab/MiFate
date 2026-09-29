# MiFate

MiFate predicts the relative response of microbial genomes to secondary treatment (SEC) and disinfection (DIS). For each stage, it reports a score and one of three labels: **Relative attenuation**, **Relative retention**, or **Uncertain**.

![MiFate workflow from microbial genomes to treatment-response predictions](assets/mifate_workflow.png)

## Run from genome FASTA files

Provide one `.fa`, `.fna`, or `.fasta` file per genome (a folder or ZIP is also accepted). File names must match the `genome_id` in a quality table:

```csv
genome_id,completeness,contamination
MAG_001,82.4,2.1
MAG_002,50.0,9.9
```

Genomes need **completeness ≥50%** and **contamination <10%**. Values are percentages.

Install Python dependencies with `python -m pip install -r requirements.txt`. Also install [Prodigal](https://github.com/hyattpd/Prodigal) and [KOfamScan](https://github.com/takaram/kofam_scan) so `prodigal` and `exec_annotation` are available. Download the KOfam profiles and their matching `ko_list` threshold file.

From the repository root, run:

```bash
bash bin/run_mifate_full.sh \
  -i /data/genomes -q /data/genome_quality.csv \
  -d /opt/kofam/profiles -k /opt/kofam/ko_list \
  -o /data/mifate_results -t 8
```

The output folder must be new or empty. The pipeline predicts proteins, accepts only KOfamScan hits marked `*` above KO-specific thresholds, builds the 8,926-feature binary KO profile, and applies the saved SEC and DIS models. **It does not retrain them.**

## Results

The main table is `mifate_results/predictions/MiFate_predictions.csv`, with genome ID, SEC/DIS scores, predictions, and applicability-domain status. The detailed table adds QC values, domain distances, and reasons for uncertain predictions. Genomes that fail QC or have invalid KO profiles appear in `excluded_genomes.csv`.

`Uncertain` means the genome is outside a model's applicability domain or its score falls between that model's class thresholds. Predictions describe **relative treatment-response tendencies**, not absolute survival in a particular treatment plant.

## If you already have a KO matrix

Provide a binary CSV with genome IDs in the first column and KO IDs in the remaining columns:

```bash
python src/mifate_predict.py \
  --ko /data/KO_matrix.csv \
  --quality /data/genome_quality.csv \
  --output /data/mifate_matrix_results
```

Missing model KOs are set to zero; extra KOs are ignored. The matrix should use the same KO-calling approach as the training data.

## Model version

This package contains the **2026-09-29 retrained models**; their predictions can differ from the earlier repository release, which remains available in the Git history. The original training matrix lacks raw KOfamScan hit records, so exact annotation compatibility for new FASTA inputs still needs verification. See [MODEL_CARD.md](MODEL_CARD.md) for model details and this limitation.
