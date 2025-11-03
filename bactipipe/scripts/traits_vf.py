# bactipipe/scripts/traits_vf.py
import os
import csv
import subprocess
from typing import List, Dict, Callable, Optional
from importlib.resources import files as resource_files

from bactipipe.scripts.core import env_cmd, set_tmp_env
from bactipipe.scripts import traits_db

# Unified output schema for per-sample VF TSV
_SCHEMA = [
    "sample","gene","allele","category","identity","coverage","length",
    "contig","start","end","strand","source","database","accession","tool"
]

# ----------------------------
# Utilities
# ----------------------------
def _log(logfn: Optional[Callable], msg: str):
    try:
        logfn(msg, "Norm")
    except Exception:
        pass

def _run(cmd: List[str], *, env: Optional[dict] = None) -> str:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr}")
    return p.stdout

def _num(x) -> float:
    try:
        return float(str(x).replace("%", "").strip())
    except Exception:
        return 0.0

# ----------------------------
# Parsers
# ----------------------------
def _vf_map() -> Dict[str,str]:
    path = resource_files("bactipipe.data").joinpath("vf_category_map.tsv")
    mp: Dict[str,str] = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = (row.get("gene") or "").strip()
            if g:
                cat = (row.get("functional_category") or "Misc").strip() or "Misc"
                mp[g] = cat
                mp[g.lower()] = cat    # case-insensitive key
    return mp

def _vf_category(gene: str, mp: Dict[str,str]) -> str:
    if not gene:
        return "Misc"
    g = gene.strip()
    # exact, lower, then base token before '-', '/', or space
    return (mp.get(g)
            or mp.get(g.lower())
            or mp.get(g.split("-",1)[0])
            or mp.get(g.lower().split("-",1)[0])
            or mp.get(g.split("/",1)[0])
            or mp.get(g.lower().split("/",1)[0])
            or mp.get(g.split(" ",1)[0])
            or mp.get(g.lower().split(" ",1)[0])
            or "Misc")

def _parse_vfinder_results(out_dir: str, sample: str, catmap: Dict[str,str]) -> List[Dict[str,str]]:
    """
    Parse VirulenceFinder extended output (results_tab.tsv produced with -x).
    Expected headers (from your example):
      Database, Virulence factor, Identity, Query / Template length,
      Contig, Position in contig, Protein function, Accession number
    """
    tsv = os.path.join(out_dir, "results_tab.tsv")
    if not os.path.exists(tsv):
        return []

    def _parse_qt(s: str) -> tuple:
        # "623 / 621" -> (623, 621)
        try:
            q, t = [x.strip() for x in s.split("/", 1)]
            return (float(q), float(t)) if (q and t) else (0.0, 0.0)
        except Exception:
            return (0.0, 0.0)

    def _parse_pos(s: str) -> tuple:
        # "1989087..1989665" -> ("1989087","1989665")
        try:
            a, b = [x.strip() for x in s.split("..", 1)]
            return (a, b)
        except Exception:
            return ("", "")

    rows: List[Dict[str,str]] = []
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            gene = (rec.get("Virulence factor") or "").strip()
            if not gene:
                continue

            ident = (rec.get("Identity") or "").strip()
            qlen, tlen = _parse_qt(rec.get("Query / Template length") or "")
            cov = f"{(qlen / tlen * 100.0):.2f}" if tlen > 0 else ""

            contig = (rec.get("Contig") or "").strip()
            start, end = _parse_pos(rec.get("Position in contig") or "")
            acc = (rec.get("Accession number") or "").strip()
            length = str(int(qlen)) if qlen > 0 else ""

            rows.append(dict(
                sample=sample,
                gene=gene,
                allele="",  # VF doesn’t report allele here
                category=_vf_category(gene, catmap),
                identity=str(ident),
                coverage=str(cov),
                length=length,
                contig=contig,
                start=start,
                end=end,
                strand="",  # not provided
                source="VirulenceFinder",
                database="virulencefinder_db",
                accession=acc,
                tool="VirulenceFinder",
            ))
    return rows


def _parse_abricate_vfdb(tsv_path: str, sample: str, catmap: Dict[str,str]) -> List[Dict[str,str]]:
    """Parse ABRicate VFDB TSV."""
    out: List[Dict[str,str]] = []
    if not os.path.exists(tsv_path):
        return out
    with open(tsv_path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for rec in rdr:
            gene = rec.get("GENE") or rec.get("PRODUCT") or rec.get("SEQUENCE", "")
            out.append(dict(
                sample=sample,
                gene=str(gene),
                allele=rec.get("PRODUCT",""),
                category=catmap.get(gene, "Misc"),
                identity=rec.get("%IDENTITY",""),
                coverage=rec.get("%COVERAGE",""),
                length=rec.get("LENGTH",""),
                contig=rec.get("SEQID",""),
                start=rec.get("START",""),
                end=rec.get("END",""),
                strand=rec.get("STRAND",""),
                source="ABRicate:vfdb",
                database="vfdb",
                accession=rec.get("ACCESSION",""),
                tool="ABRicate",
            ))
    return out


# ----------------------------
# Main entry
# ----------------------------
def run_vf_for_sample(
    sample: str,
    fasta: str,
    out_root: str,
    *,
    min_id: float,
    min_cov: float,
    enable_vfinder: bool,
    enable_vfdb: bool,
    threads: int,
    logger_fn: Optional[Callable] = None,
    db_root: Optional[str] = None,   # <-- NEW: accept DB root passed from run_traits
) -> List[Dict[str,str]]:
    """
    Execute virulence detection for one sample.
    - DTU VirulenceFinder (Python tool) if enable_vfinder=True
    - ABRicate/VFDB if enable_vfdb=True
    Preference: keep VirulenceFinder hits when overlapping with VFDB.
    """
    catmap = _vf_map()
    rows: List[Dict[str,str]] = []

    raw_root = os.path.join(out_root, "raw")
    os.makedirs(raw_root, exist_ok=True)

    # Per-sample temp directory so all tools write into out_root/.tmp/...
    tmp_dir = os.path.join(out_root, ".tmp", "vf", sample)
    os.makedirs(tmp_dir, exist_ok=True)
    env = set_tmp_env(tmp_dir)

    # ---------------- VirulenceFinder (DTU) ----------------
    if enable_vfinder:
        vfinder_out = os.path.join(raw_root, "virulencefinder", sample)
        os.makedirs(vfinder_out, exist_ok=True)
        sentinel = os.path.join(vfinder_out, ".done")
        if not os.path.exists(sentinel):
            _log(logger_fn, f"[{sample}] Running VirulenceFinder (DTU Python tool)...")
            # Resolve DB path from unified DB root (BACTIPIPE_DB_DIR or --db-dir)
            vf_db = traits_db.find_virulencefinder_db(db_root)
            cmd = env_cmd("virulencefinder") + ["-i", fasta, "-o", vfinder_out, "-x"]
            if vf_db:
                cmd.extend(["-p", vf_db])  # modern virulencefinder supports -p for db path
            _run(cmd, env=env)
            open(sentinel, "w").close()
        rows += _parse_vfinder_results(vfinder_out, sample, catmap)

    # ---------------- ABRicate/VFDB ----------------
    if enable_vfdb:
        vfdb_dir = os.path.join(raw_root, "abricate", "vfdb")
        os.makedirs(vfdb_dir, exist_ok=True)
        vfdb_tsv = os.path.join(vfdb_dir, f"{sample}.tsv")
        if not os.path.exists(vfdb_tsv):
            _log(logger_fn, f"[{sample}] Running abricate --db vfdb ...")
            cmd = env_cmd("abricate") + ["--db", "vfdb", fasta]
            tsv = _run(cmd, env=env)
            with open(vfdb_tsv, "w") as f:
                f.write(tsv)
        rows += _parse_abricate_vfdb(vfdb_tsv, sample, catmap)

    # ---------------- Merge & preference ----------------
    # Prefer VirulenceFinder over VFDB when hits overlap by (gene, contig, start, end)
    merged: List[Dict[str,str]] = []
    seen: Dict[tuple, int] = {}
    for r in rows:
        key = (r.get("gene",""), r.get("contig",""), r.get("start",""), r.get("end",""))
        if key in seen:
            i = seen[key]
            if merged[i]["source"] != "VirulenceFinder" and r["source"] == "VirulenceFinder":
                merged[i] = r
        else:
            seen[key] = len(merged)
            merged.append(r)

    # ---------------- Threshold filter ----------------
    kept = [r for r in merged if _num(r.get("identity","0")) >= min_id and _num(r.get("coverage","0")) >= min_cov]

    # ---------------- Write sample VF TSV ----------------
    out_tsv = os.path.join(out_root, f"{sample}.vf.tsv")
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_SCHEMA, delimiter="\t")
        w.writeheader()
        for r in kept:
            w.writerow({k: r.get(k, "") for k in _SCHEMA})

    return kept
