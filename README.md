[readme_Mifate.md](https://github.com/user-attachments/files/27628778/readme_Mifate.md)
# MiFate: Microbial Fate Predictor for Wastewater Treatment

MiFate is a genome-informed pipeline for predicting the persistence of microbial genomes (MAGs) in wastewater treatment plants (WWTPs). It integrates protein prediction, KEGG annotation, KO matrix construction, and model-based prediction to determine which MAGs are likely removable or persistent.

## Repository Structure

```
MiFate/
├─ bin/                    # Executable scripts
│  ├─ run_mifate_full.sh        # One-command full pipeline
│  ├─ run_prodigal_kofam.sh    # Protein prediction & KEGG annotation
│  ├─ build_binary_ko_matrix.py # Convert KOfam mapper to binary matrix
│  └─ predict.py               # XGBoost prediction script
├─ models/                  # Pre-trained models and feature files
│  ├─ SEC_model.xgb
│  ├─ DIS_model.xgb
│  ├─ SEC_features.txt
│  ├─ DIS_features.txt
│  └─ Best_thresholds.tsv
├─ db/                      # KOfam database
├─ input/                   # Example MAG input folder
├─ output/                  # Output folder
└─ README.md
```

## Prerequisites

- Linux or macOS with bash
- [Prodigal](https://bioconda.github.io/recipes/prodigal/README.html) for protein prediction
- [KOfamScan](https://github.com/takaram/kofam_scan) for KEGG annotation (`exec_annotation` + `profiles/` + `ko_list`)
- Python 3.8+ with the following packages:
```bash
pip install pandas numpy xgboost
```

## Running MiFate

Use the main script to run the full pipeline:
```bash
bash bin/run_mifate_full.sh -i /path/to/MAGs -o /path/to/output -d /path/to/kofam_db --delta 0.05 --min_active 3
```

**Parameters:**
- `-i`: Input MAG directory (each MAG as `.fna`)
- `-o`: Output directory
- `-d`: KOfamScan database directory
- `--delta`: Gray-zone half-width around prediction threshold (default 0.05)
- `--min_active`: Minimum active features per MAG (default 3)
- `-h, --help`: Show help message

## Outputs

- `MAG_KEGG_binary_matrix.csv`: Binary matrix of MAGs vs. KEGG orthologs
- `MAG_pred_SEC_DIS.csv`: Prediction results with labels (`hard`, `easy`, `uncertain`)
- `quadrant_counts_with_uncertain.tsv`, `SEC_class_counts.tsv`, `DIS_class_counts.tsv`: Summary statistics of predictions

## Notes

- Predictions are per-MAG independent; results do not depend on the number of input MAGs.
- Start with a small set of MAGs for testing.
- Ensure KOfamScan database and Prodigal are installed and accessible.
- The repository provides `bin/` and `models/` folders; users only need to provide their MAGs and KOfam database.

## Reference

If using this pipeline, please cite: [xxxx].
