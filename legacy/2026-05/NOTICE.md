# Previous MiFate GitHub release

These are the scripts, model files, and README from the repository's previous
`main` revision (commit `60cf58df6a4598057f5511a9b681f79cc5fd2d21`). They
are retained to make historical runs traceable. The current deployment is at
the repository root. Do not combine the old model weights or feature lists
with the new score thresholds or applicability-domain reference data.

The old KOfamScan call used `-f mapper` with `profiles/` and `ko_list` supplied
separately. KOfamScan's mapper format reports above-threshold KO assignments,
so the old annotation step already applied the official KO threshold rule.
