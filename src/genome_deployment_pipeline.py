"""MiFate genome-to-prediction deployment pipeline.

One FASTA file represents one genome (it may contain many contigs).  The
pipeline runs Prodigal, KOfamScan, builds the frozen 8,926-column KO matrix,
and delegates SEC/DIS scoring to mifate_predict.py.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, shutil, subprocess, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(__file__).resolve().parent

def run(cmd: list[str], log: Path, cwd: Path | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        p = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd[0]}. See {log}")

def read_fasta_ids(path: Path) -> list[str]:
    ids=[]
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                x=line[1:].strip().split()[0]
                if x: ids.append(x)
    if not ids: raise ValueError(f"No FASTA records found: {path}")
    return ids

def genome_files(inp: Path, work: Path) -> list[Path]:
    if inp.is_file() and inp.suffix.lower()==".zip":
        dst=work/"genomes"; dst.mkdir()
        with zipfile.ZipFile(inp) as z:
            for entry in z.infolist():
                p=Path(entry.filename)
                if p.suffix.lower() in {".fna",".fa",".fasta"} and not entry.is_dir():
                    out=dst/p.name
                    if out.exists():
                        raise ValueError(f"Duplicate FASTA filename inside ZIP: {p.name}")
                    with z.open(entry) as source, out.open("wb") as target:
                        shutil.copyfileobj(source,target)
        fs=sorted(dst.glob("*"))
    elif inp.is_dir(): fs=sorted([p for p in inp.iterdir() if p.suffix.lower() in {".fna",".fa",".fasta"}])
    elif inp.is_file(): fs=[inp]
    else: raise FileNotFoundError(inp)
    if not fs: raise ValueError("No .fna/.fa/.fasta files found")
    for f in fs: read_fasta_ids(f)
    return fs

def read_quality(path: Path) -> dict[str, tuple[float,float]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows=list(csv.DictReader(fh, delimiter="\t" if path.suffix.lower() in {".tsv",".txt"} else ","))
    if not rows: raise ValueError("Quality table is empty")
    def col(row, names):
        for k in row:
            if k.strip().lower() in names: return row[k]
        return None
    out={}
    for r in rows:
        gid=col(r,{"genome_id","genome","mag","mag_id","id"})
        c=col(r,{"completeness","completeness_percent"}); x=col(r,{"contamination","contamination_percent"})
        if not gid or c is None or x is None: raise ValueError("Quality table requires genome_id, completeness, contamination")
        if gid in out: raise ValueError(f"Duplicate genome_id in quality table: {gid}")
        try:
            out[gid]=(float(c),float(x))
        except ValueError as exc:
            raise ValueError(f"Non-numeric quality value for {gid}") from exc
    return out

def parse_kofamscan(path: Path) -> set[str]:
    """Accept only KOfamScan detail-tsv hits marked '*' above its KO threshold."""
    kos=set()
    starred=re.compile(r"^\s*\*\s+(\S+)\s+(K\d{5})(?:\s|$)")
    header_found=False
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip(): continue
        if line.lstrip().startswith("#"):
            if "gene name" in line and "KO" in line and "score" in line:
                header_found=True
            continue
        if not header_found:
            raise ValueError(f"Expected KOfamScan detail-tsv header before hits: {path}:{lineno}")
        if line.lstrip().startswith("*"):
            m=starred.match(line)
            if not m:
                raise ValueError(f"Malformed starred KOfamScan hit: {path}:{lineno}")
            kos.add(m.group(2))
    if not header_found:
        raise ValueError(f"Expected KOfamScan detail-tsv output with header: {path}")
    return kos

def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", required=True, type=Path, help="One .fna, a directory of genomes, or a .zip")
    ap.add_argument("--quality", required=True, type=Path, help="CSV/TSV with genome_id, completeness, contamination")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--prodigal", default="prodigal")
    ap.add_argument("--kofamscan", default="exec_annotation")
    ap.add_argument("--kofam-profile", type=Path, required=True,
                    help="KOfam profile HMM directory, .hmm file, or .hal file")
    ap.add_argument("--kofam-ko-list", type=Path, required=True,
                    help="KOfam ko_list with official KO-specific score thresholds")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--artifacts", type=Path, default=ROOT/"models/applicability_domain")
    ap.add_argument("--rules", type=Path, default=ROOT/"models/deployment_rules.json")
    ap.add_argument("--keep-intermediates", action="store_true")
    a=ap.parse_args()
    if a.threads < 1:
        ap.error("--threads must be positive")
    if not a.kofam_profile.exists() or not (a.kofam_profile.is_dir() or a.kofam_profile.suffix.lower() in {".hmm", ".hal"}):
        ap.error("--kofam-profile must be a KOfam profiles directory, .hmm, or .hal file")
    if not a.kofam_ko_list.is_file():
        ap.error("--kofam-ko-list must be the KOfam ko_list file")
    out=a.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        ap.error("--output must be a new or empty directory to avoid mixing runs")
    out.mkdir(parents=True, exist_ok=True)
    work=out/"work"; work.mkdir(exist_ok=True); log=out/"pipeline.log"
    q=read_quality(a.quality); fs=genome_files(a.genomes,work)
    ids=[p.stem for p in fs]
    if len(ids)!=len(set(ids)): raise ValueError("Genome filenames must be unique after removing extensions")
    missing=[g for g in ids if g not in q]
    if missing: raise ValueError("Quality table lacks genome_id: "+", ".join(missing))
    model_cols=[]
    for s in ("SEC","DIS"):
        ref=json.loads((a.artifacts/s/"reference.json").read_text(encoding="utf-8"))
        model_cols.append(ref["n_features"])
    if len(set(model_cols))!=1 or model_cols[0]!=8926: raise ValueError(f"Frozen model feature count unexpected: {model_cols}")
    import numpy as np
    with np.load(a.artifacts/"SEC"/"reference.npz",allow_pickle=False) as z:
        cols=[str(x) for x in z["features"].tolist()]
    with np.load(a.artifacts/"DIS"/"reference.npz",allow_pickle=False) as z:
        dis_cols=[str(x) for x in z["features"].tolist()]
    if cols != dis_cols:
        raise ValueError("SEC and DIS feature orders differ; cannot build one shared KO matrix")
    kos={}
    qc=[]
    for f in fs:
        gid=f.stem; comp,cont=q[gid]
        passed=(math.isfinite(comp) and math.isfinite(cont)
                and 50<=comp<=100 and 0<=cont<10)
        qc.append({"genome_id":gid,"completeness":comp,"contamination":cont,"QC_pass":passed})
        if not passed: continue
        gdir=work/gid; gdir.mkdir(parents=True,exist_ok=True); faa=gdir/(gid+".faa"); kof=gdir/(gid+".kofam.tsv")
        run([a.prodigal,"-i",str(f),"-a",str(faa),"-p","meta","-q"],log)
        run([a.kofamscan,"--profile",str(a.kofam_profile),"--ko-list",str(a.kofam_ko_list),
             "--format","detail-tsv","--tmp-dir",str(gdir/"kofam_tmp"),
             "-o",str(kof),"--cpu",str(a.threads),str(faa)],log)
        kos[gid]=parse_kofamscan(kof)
    # Build exact model-column matrix from the persisted feature order.
    # reference.npz is authoritative for order; avoid trusting a separate input matrix.
    matrix=work/"MiFate_KO_matrix.csv"
    with matrix.open("w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh); w.writerow(["MAG_ID"]+cols)
        for gid in ids: w.writerow([gid]+[1 if k in kos.get(gid,set()) else 0 for k in cols])
    quality=work/"MiFate_quality.csv"
    with quality.open("w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh); w.writerow(["MAG_ID","Completeness","Contamination"])
        for r in qc: w.writerow([r["genome_id"],r["completeness"],r["contamination"]])
    predictor=SRC/"mifate_predict.py"
    py=sys.executable
    run([py,str(predictor),"--ko",str(matrix),"--quality",str(quality),"--output",str(out/"predictions"),"--artifacts",str(a.artifacts),"--rules",str(a.rules)],log)
    (out/"qc_summary.csv").write_text("genome_id,completeness,contamination,QC_pass\n"+"\n".join(f'{r["genome_id"]},{r["completeness"]},{r["contamination"]},{r["QC_pass"]}' for r in qc)+"\n",encoding="utf-8")
    with (out/"predictions"/"MiFate_predictions.csv").open(newline="",encoding="utf-8-sig") as fh:
        n_predicted=sum(1 for _ in csv.DictReader(fh))
    def sha256(path):
        digest=hashlib.sha256()
        with Path(path).open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024*1024),b""):
                digest.update(chunk)
        return digest.hexdigest()
    manifest={"model_release":"MiFate retrained 2026-09-29 (not the original manuscript model)","input":str(a.genomes.resolve()),"quality":str(a.quality.resolve()),"n_genomes":len(ids),"n_qc_pass":sum(r["QC_pass"] for r in qc),"n_predicted":n_predicted,"n_excluded":len(ids)-n_predicted,"feature_count":len(cols),"kofam_hit_policy":"Only '*' marked detail-tsv hits at KO-specific thresholds","kofam_profile":str(a.kofam_profile.resolve()),"kofam_ko_list":str(a.kofam_ko_list.resolve()),"kofam_ko_list_sha256":sha256(a.kofam_ko_list),"n_accepted_kos_by_genome":{gid:len(ks) for gid,ks in kos.items()},"outputs":["predictions/MiFate_predictions.csv","predictions/MiFate_predictions_detailed.csv","predictions/excluded_genomes.csv","predictions/run_audit.json","qc_summary.csv","pipeline.log"]}
    (out/"run_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    if not a.keep_intermediates: shutil.rmtree(work,ignore_errors=True)
    print(json.dumps(manifest,indent=2))

if __name__=="__main__": main()


