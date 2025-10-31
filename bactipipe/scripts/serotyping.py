# bactipipe/scripts/serotyping.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, Tuple, List
import csv
import logging
import shutil

from .runner import run as _run
from .runner import ensure_dir as _ensure_dir

def _ensure_prefix(val: str, prefix: str) -> str:
    v = (val or "").strip()
    if not v or v.lower() == "unknown":
        return "Unknown"
    up = v.upper()
    return v if up.startswith(prefix.upper()) else f"{prefix}{v}"

def run_serotypefinder(
    genome: Path,
    outdir: Path,
    serodb_dir: Path,
    logger: logging.Logger,
    env_name: Optional[str] = None,
) -> Tuple[str, str]:
    """
    Run serotypefinder.py and return (O_type, H_type) if present.
    If multiple O hits, choose the one with full-length coverage (>=1.0) and highest identity.
    Returns 'Unknown' when not found.
    """
    _ensure_dir(outdir)
    cmd = ["serotypefinder.py", "-i", str(genome), "-o", str(outdir), "-p", str(serodb_dir), "-x"]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (outdir / "serotypefinder.stdout.txt").write_text(out)
    if rc != 0:
        logger.warning(f"serotypefinder failed for {genome.name}")
        return "Unknown", "Unknown"

    tsv = outdir / "results_tab.tsv"
    if not tsv.exists():
        logger.warning(f"serotypefinder results_tab.tsv missing for {genome.name}")
        return "Unknown", "Unknown"

    O_type = "Unknown"
    H_type = "Unknown"
    best_O_ident = -1.0

    try:
        for row in tsv.read_text().splitlines():
            if not row.strip():
                continue
            cols = row.split("\t")
            if cols[0] == "Database":
                continue
            db = cols[0]
            if db == "H_type":
                if len(cols) >= 3 and cols[2].strip():
                    H_type = cols[2].strip()
            elif db == "O_type" and len(cols) >= 5:
                try:
                    ident = float(cols[3])
                    cov = cols[4].split(" / ")
                    ratio = (int(cov[0]) / int(cov[1])) if len(cov) == 2 and int(cov[1]) else 0.0
                except Exception:
                    ident, ratio = 0.0, 0.0
                this_O = cols[2].strip() if len(cols) >= 3 else ""
                if this_O and ratio >= 1.0 and ident >= best_O_ident:
                    best_O_ident = ident
                    O_type = this_O
        return O_type or "Unknown", H_type or "Unknown"
    except Exception as e:
        logger.warning(f"Failed to parse serotypefinder output for {genome.name}: {e}")
        return "Unknown", "Unknown"

def run_seqsero2(
    genome: Path,
    outdir: Path,
    logger: logging.Logger,
    threads: int = 4,
    env_name: Optional[str] = None,
) -> Tuple[str, str, str]:
    """
    Return (serotype, antigenic_formula, species). On failure, ('NA','NA','NA').
    Tries k-mer mode first; if the known get_input_K/input_Ks bug appears, falls back to assembly mode.
    Threads are clamped to 1–5 to satisfy SeqSero2 limits.
    """
    _ensure_dir(outdir)

    try:
        t = max(1, min(5, int(threads)))
    except Exception:
        t = 4

    def _call(mode: str) -> Tuple[int, str]:
        cmd = ["SeqSero2_package.py", "-m", mode, "-t", str(t), "-i", str(genome), "-d", str(outdir)]
        rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
        with (outdir / "seqsero2.stdout.txt").open("a") as fh:
            fh.write(f"\n=== MODE {mode} (t={t}) ===\n")
            fh.write(out)
        return rc, out

    rc, out = _call("k")
    if rc != 0 and ("get_input_K" in out or "local variable 'input_Ks' referenced before assignment" in out):
        logger.info("SeqSero2 k-mer mode hit get_input_K bug; retrying with assembly mode (-m a).")
        rc, out = _call("a")

    if rc != 0:
        logger.warning(f"SeqSero2 failed for {genome.name}; see {outdir/'seqsero2.stdout.txt'}")
        return "NA", "NA", "NA"

    candidates = [outdir / "SeqSero_result.tsv", outdir / "results_tab.tsv"]
    tsv = next((p for p in candidates if p.exists()), None)
    if tsv is None:
        logger.warning(f"SeqSero2 output TSV not found for {genome.name}")
        return "NA", "NA", "NA"

    try:
        rows = [line.split("\t") for line in tsv.read_text().splitlines() if line.strip()]
        if not rows:
            return "NA", "NA", "NA"
        header = rows[0]
        data = rows[1] if len(rows) > 1 else None
        if not data:
            return "NA", "NA", "NA"
        colmap = {h.lower(): i for i, h in enumerate(header)}
        serotype = data[colmap.get("predicted serotype", -1)] if colmap.get("predicted serotype", -1) != -1 else "NA"
        formula  = data[colmap.get("predicted antigenic profile", -1)] if colmap.get("predicted antigenic profile", -1) != -1 else "NA"
        species  = data[colmap.get("species", -1)] if colmap.get("species", -1) != -1 else "NA"
        if serotype == "NA" and len(data) >= 10:
            serotype = data[8] or serotype
        if formula == "NA" and len(data) >= 8:
            formula = data[7]
        return serotype or "NA", formula or "NA", species or "NA"
    except Exception as e:
        logger.warning(f"Failed to parse SeqSero2 TSV for {genome.name}: {e}")
        return "NA", "NA", "NA"

def run_kleborate(
    genome: Path,
    outdir: Path,
    logger: logging.Logger,
    env_name: Optional[str] = None,
) -> Tuple[str, str, str]:
    """
    Run Kleborate on a single assembly and parse K and O locus/type.
    Returns (K_call, O_call, species_label). On failure -> ('Unknown','Unknown','Klebsiella').
    """
    outdir.mkdir(parents=True, exist_ok=True)
    if not shutil.which("kleborate"):
        logger.warning("kleborate not found on PATH; falling back to Unknown.")
        return "Unknown", "Unknown", "Klebsiella"

    cmd = ["kleborate", "-a", str(genome), "-o", str(outdir), "-p", "kpsc", "--trim_headers"]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (outdir / "kleborate.stdout.txt").write_text(out)
    if rc != 0:
        logger.warning(f"Kleborate failed for {genome.name}")
        return "Unknown", "Unknown", "Klebsiella"

    candidates = sorted(outdir.glob("*output.txt"))
    if not candidates:
        logger.warning(f"Kleborate output file not found for {genome.name}")
        return "Unknown", "Unknown", "Klebsiella"
    tsv = next((p for p in candidates if "kleb" in p.name.lower()), candidates[0])

    try:
        rows = [r for r in csv.reader(tsv.open("r"), delimiter="\t") if any(c.strip() for c in r)]
        if not rows:
            return "Unknown", "Unknown", "Klebsiella"

        header = [h.strip() for h in rows[0]]
        hmap = {h.strip().lower(): i for i, h in enumerate(header)}
        data = rows[1] if len(rows) > 1 else None
        if not data:
            return "Unknown", "Unknown", "Klebsiella"

        def get_any(keys: List[str]) -> str:
            for k in keys:
                if k in hmap and hmap[k] < len(data):
                    val = (data[hmap[k]] or "").strip()
                    if val:
                        return val
            return ""

        species = get_any(["species", "klebsiella_pneumo_complex__species", "kpsc__species"]) or "Klebsiella"
        K_call  = get_any(["k_locus", "k_type", "kpsc__k_locus", "kpsc__k_type"])
        O_call  = get_any(["o_locus", "o_type", "kpsc__o_locus", "kpsc__o_type"])
        # log suboptimal confidence if present
        for lab, conf in (("K", get_any(["k_locus_confidence", "kpsc__k_locus_confidence"])),
                          ("O", get_any(["o_locus_confidence", "kpsc__o_locus_confidence"]))):
            if conf and conf.lower() not in {"very high", "high", "good", "typeable"}:
                logger.info(f"Kleborate {lab} confidence ({conf}) for {genome.name}")

        return (K_call or "Unknown"), (O_call or "Unknown"), (species or "Klebsiella")
    except Exception as e:
        logger.warning(f"Failed to parse Kleborate output for {genome.name}: {e}")
        return "Unknown", "Unknown", "Klebsiella"

def serotype_dispatch(
    organism: str,
    genome: Path,
    outdir: Path,
    db_root: Path,
    logger: logging.Logger,
    threads: int = 4,
    envs: Optional[Dict[str, Optional[str]]] = None,
) -> Tuple[str, str, str]:
    """
    Return (serotype, formula, species).
    - Salmonella: SeqSero2 -> (serotype, antigenic_formula, species)
    - E. coli:    SerotypeFinder -> (Oxx:Hyy, Oxx:Hyy, 'Escherichia coli')
    - Klebsiella: Kleborate; fallback SerotypeFinder O-only; serotype==formula
    - Default:    SerotypeFinder -> (Oxx:Hyy, Oxx:Hyy, <organism>)
    """
    envs = envs or {}
    org = organism.lower()

    if org == "salmonella":
        sdir = outdir / "seqsero2"
        sero, formula, species = run_seqsero2(genome, sdir, logger, threads=threads, env_name=envs.get("seqsero2"))
        return sero, formula, species if species != "NA" else "Salmonella"

    if org in {"ecoli", "e. coli", "escherichia", "escherichia coli"}:
        sdir = outdir / "serotypefinder"
        serodb = db_root / "serotypefinder_db"
        O_type, H_type = run_serotypefinder(genome, sdir, serodb, logger, env_name=envs.get("serotypefinder"))
        o = _ensure_prefix(O_type, "O")
        h = _ensure_prefix(H_type, "H")
        combined = "Unknown" if (o == "Unknown" and h == "Unknown") else (f"{o}:{h}" if (o != "Unknown" and h != "Unknown") else (o if o != "Unknown" else h))
        return combined, combined, "Escherichia coli"

    if org.startswith("kleb"):
        kdir = outdir / "kleborate"
        K_call, O_call, species = run_kleborate(genome, kdir, logger, env_name=envs.get("kleborate"))
        if K_call == "Unknown" and O_call == "Unknown":
            sdir = outdir / "serotypefinder"
            serodb = db_root / "serotypefinder_db"
            O_type, _ = run_serotypefinder(genome, sdir, serodb, logger, env_name=envs.get("serotypefinder"))
            o = _ensure_prefix(O_type, "O")
            combined = o if o != "Unknown" else "Unknown"
            return combined, combined, species or "Klebsiella"
        parts = []
        if K_call and K_call != "Unknown": parts.append(K_call)  # e.g., KL2
        if O_call and O_call != "Unknown": parts.append(O_call)  # e.g., O1
        combined = ":".join(parts) if parts else "Unknown"
        return combined, combined, species or "Klebsiella"

    # default
    sdir = outdir / "serotypefinder"
    serodb = db_root / "serotypefinder_db"
    O_type, H_type = run_serotypefinder(genome, sdir, serodb, logger, env_name=envs.get("serotypefinder"))
    o = _ensure_prefix(O_type, "O")
    h = _ensure_prefix(H_type, "H")
    combined = "Unknown" if (o == "Unknown" and h == "Unknown") else (f"{o}:{h}" if (o != "Unknown" and h != "Unknown") else (o if o != "Unknown" else h))
    return combined, combined, organism
