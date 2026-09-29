# Model card: MiFate retrained 2026-09-29

This repository contains the **retrained** MiFate SEC and DIS XGBoost
classifiers. It is a reproducible new model version, not the manuscript's
earlier frozen model. Its model weights, Train80 reference genomes, domain
cutoffs, and score thresholds differ from the earlier version.

| Item | Value |
| --- | --- |
| Features | 8,926 binary KOs, in the saved feature order |
| Training source | 7,890 MAG KO profiles and wastewater treatment response data |
| Output | SEC and DIS score; relative attenuation, relative retention, or uncertain |
| Applicability domain | Mean binary Jaccard distance to five nearest Train80 genomes, stage-specific 99th percentile cutoff |
| SEC cutoffs | attenuation ≤0.4990727305; retention ≥0.5009139776 |
| DIS cutoffs | attenuation ≤0.4991326034; retention ≥0.8531948328 |
| Held-out AUROC | SEC 0.9683; DIS 0.9038 |

The threshold selection used Train80 out-of-fold predictions and an empirical
accepted-class error ceiling of 25%, with at least 30 accepted genomes per
class. `validation/training_summary.json` and
`validation/family_bootstrap_95CI.json` contain the source summaries. The
model reference bundles verify checksums of `model.json` and `reference.npz`
at prediction time.

**Annotation compatibility remains an open release check.** The binary
training and earlier WHO KO matrices lack original KOfamScan scores and `*`
markers; their KO calling provenance cannot be established from these matrices
alone. This repository's FASTA pipeline only accepts official starred hits.
Before interpreting new-genome predictions against this model, verify that
training KO calls were generated with comparable Prodigal and KOfamScan
settings and KOfam database release. If the training matrix included
below-threshold calls, regenerate training features and retrain. Reprocessing
an external matrix only would not repair a training mismatch.

The scores are classifier outputs for community-normalized treatment-response
labels; they are not calibrated survival probabilities. Predictions labeled
uncertain are abstentions due to domain distance or a score between cutoffs.
