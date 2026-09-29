"""Check the KO threshold boundary and a full run with the bundled models.

The FASTA integration test uses fake Prodigal and KOfamScan executables so CI
does not need the multi-gigabyte KOfam database; scoring uses the real SEC/DIS
model artifacts. A real-data annotation run is still required for release QA.
"""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from genome_deployment_pipeline import parse_kofamscan  # noqa: E402


HEADER = "# gene name\tKO\tthrshld\tscore\tE-value\tKO definition\n"


class KOfamThresholdTests(unittest.TestCase):
    def test_only_starred_hits_count(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hits.tsv"
            path.write_text(
                HEADER
                + "*\tgene1\tK00003\t20\t35\t1e-10\taccepted\n"
                + "gene2\tK00004\t20\t10\t1e-2\tbelow threshold\n"
                + "  *  gene3  K00005  20  27  1e-7  accepted\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_kofamscan(path), {"K00003", "K00005"})

    def test_wrong_format_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hits.tsv"
            path.write_text("gene1\tK00003\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_kofamscan(path)
            path.write_text(HEADER + "*\tgene1\tNO_KO\t20\t30\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_kofamscan(path)

    def test_fasta_to_real_model_with_mock_annotations(self):
        with tempfile.TemporaryDirectory() as d:
            temp = Path(d)
            genomes = temp / "genomes"
            genomes.mkdir()
            for name in ("qc_boundary", "bad_qc"):
                (genomes / f"{name}.fa").write_text(">contig\nATGAAATAA\n", encoding="utf-8")
            quality = temp / "quality.csv"
            quality.write_text(
                "genome_id,completeness,contamination\n"
                "qc_boundary,50.0,9.9\nbad_qc,49.9,1.0\n",
                encoding="utf-8",
            )
            profile = temp / "profiles"
            profile.mkdir()
            ko_list = temp / "ko_list"
            ko_list.write_text("example for fake executable\n", encoding="utf-8")
            prodigal = temp / "fake_prodigal"
            prodigal.write_text(
                "#!/usr/bin/env python3\n"
                "import sys, pathlib\n"
                "pathlib.Path(sys.argv[sys.argv.index('-a')+1]).write_text('>gene1\\nMKK\\n')\n",
                encoding="utf-8",
            )
            kofam = temp / "fake_kofam"
            kofam.write_text(
                "#!/usr/bin/env python3\n"
                "import sys, pathlib\n"
                "args=sys.argv\n"
                "assert args[args.index('--format')+1]=='detail-tsv'\n"
                "assert '--profile' in args and '--ko-list' in args\n"
                "pathlib.Path(args[args.index('-o')+1]).write_text("
                + repr(HEADER + "*\tgene1\tK00003\t20\t35\t1e-10\taccepted\n"
                       + "gene2\tK00004\t20\t10\t1e-2\tbelow threshold\n") + ")\n",
                encoding="utf-8",
            )
            prodigal.chmod(0o755)
            kofam.chmod(0o755)
            out = temp / "out"
            subprocess.run([
                sys.executable, str(ROOT / "src/genome_deployment_pipeline.py"),
                "--genomes", str(genomes), "--quality", str(quality),
                "--kofam-profile", str(profile), "--kofam-ko-list", str(ko_list),
                "--prodigal", str(prodigal), "--kofamscan", str(kofam),
                "--output", str(out), "--keep-intermediates",
            ], check=True, cwd=ROOT, capture_output=True, text=True)
            with (out / "work/MiFate_KO_matrix.csv").open(newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 2)
            good = next(row for row in rows if row["MAG_ID"] == "qc_boundary")
            self.assertEqual((good["K00003"], good["K00004"]), ("1", "0"))
            with (out / "predictions/MiFate_predictions.csv").open(
                    newline="", encoding="utf-8-sig") as fh:
                predictions = list(csv.DictReader(fh))
            self.assertEqual([r["Genome_ID"] for r in predictions], ["qc_boundary"])
            audit = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(audit["n_accepted_kos_by_genome"], {"qc_boundary": 1})
            self.assertEqual((audit["n_qc_pass"], audit["n_predicted"]), (1, 1))


if __name__ == "__main__":
    unittest.main()
