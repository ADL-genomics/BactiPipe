#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Dict
import sys
import logging
import textwrap
import subprocess
import shutil
from venv import logger
from bactipipe.scripts.caps import CAPS, supported_rows
from bactipipe.scripts import ska as ska2
from bactipipe.scripts.runner import run as _run, ensure_dir as _ensure_dir
from bactipipe.scripts.serotyping import serotype_dispatch
from bactipipe.scripts import caps
import time
from contextlib import contextmanager
from bactipipe.scripts import utils as U
# ----------------------------
# Data structures
# ----------------------------
@dataclass
class SampleRecord:
    sample: str
    isolate: str
    specimen: str
    genome_path: Path
    exists: bool
    notes: str = ""
# ----------------------------
# Utilities
# ----------------------------
@contextmanager
def _step(logger: logging.Logger, label: str):
    logger.info(f"{label} ...")
    t0 = time.time()
    try:
        yield
    finally:
        dt = time.time() - t0
        logger.info(f"{label} ✓ ({dt:.1f}s)")

def _init_logging(outdir: Path, verbose: bool = False) -> logging.Logger:
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "run.log"

    # line-buffer stdout so messages appear immediately
    try:
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass

    logger = logging.getLogger("type_genomes")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))

    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))

    logger.addHandler(ch)
    logger.addHandler(fh)

    logger.info("====== Typing run: Phase 2 (core typing) ======")
    logger.info(f"Log file: {log_file}")
    return logger


def _default_sample_sheet(accession: str) -> Optional[Path]:
    candidates: List[Path] = []
    env_dir = os.environ.get("BACTIPIPE_SAMPLE_SHEETS_DIR")
    if env_dir:
        candidates.append(Path(env_dir) / f"{accession}.txt")

    home = Path.home()
    candidates.extend(
        [
            home / "Dropbox" / "typing_samples_sheets" / f"{accession}.txt",
            home / "Dropbox" / "Genomics ADL" / "NCBI_samples" / "typing_samples_sheets" / f"{accession}.txt",
        ]
    )
    for p in candidates:
        if p.exists():
            return p
    return None


def _ensure_layout(outdir: Path, create_tmp: bool = True) -> None:
    subdirs = [
        outdir / "serotype",
        outdir / "mlst",
        outdir / "cgMLST",
        outdir / "ani",
        outdir / "skani",         # <- add this
    ]
    if create_tmp:
        subdirs.append(outdir / "tmp")
    for d in subdirs:
        d.mkdir(parents=True, exist_ok=True)

def _check_tool(
    name: str,
    logger: logging.Logger,
    env_name: Optional[str] = None,
    test_args: Optional[Tuple[str, ...]] = None,
) -> bool:
    """
    True if the tool can be invoked. When env_name is provided, we test it via conda/mamba run
    (no shutil.which pre-check). Many tools return non-zero for -h, so we only fail on rc==127
    or clear 'not found' messages.
    """
    default_tests: Dict[str, Tuple[str, ...]] = {
        "mlst.py": ("-h",),
        "cgMLST.py": ("-h",),
        "serotypefinder.py": ("-h",),
        "SeqSero2_package.py": ("-h",),
        "kleborate": ("--version",),
        "skani": ("--version",),
        "fastANI": ("--version",),
        "ska": ("--version",),  # SKA2
    }
    args = test_args or default_tests.get(name, ("--version",))

    # If an env is requested, probe inside that env (do NOT rely on current PATH)
    if env_name:
        rc, out = _run([name] + list(args), cwd=None, log=logger, env_name=env_name)
        txt = (out or "").lower()
        if rc == 127 or "not found" in txt or "command not found" in txt or "no such file" in txt:
            logger.error(f"Required tool not found in env '{env_name}': {name}")
            return False
        return True

    # No env: require it on current PATH, then try a lightweight probe
    if shutil.which(name) is None:
        logger.error(f"Required tool not found on PATH: {name}")
        return False
    rc, out = _run([name] + list(args), cwd=None, log=logger, env_name=None)
    if rc == 127:
        logger.error(f"{name}: execution failed (command not found).")
        return False
    return True

def _read_sample_sheet(path: Path, logger: logging.Logger) -> List[Tuple[str, str, str]]:
    rows: List[Tuple[str, str, str]] = []
    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {path}")
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header_checked = False
        for raw in reader:
            if not raw or all((c or "").strip() == "" for c in raw):
                continue
            cells = [c.strip() for c in raw]
            if not header_checked:
                header_checked = True
                lower = [c.lower() for c in cells[:4]]
                # accept either 3-col or legacy 4-col header
                if lower[:3] == ["sample", "isolate", "specimen"]:
                    logger.debug("Detected header row in sample sheet; skipping it.")
                    continue
            if len(cells) < 3:
                logger.warning(f"Skipping malformed line (expected ≥3 columns): {raw}")
                continue
            sample, isolate, specimen = cells[:3]
            rows.append((sample, isolate, specimen))
    if not rows:
        raise ValueError("Sample sheet appears empty after parsing valid rows.")
    logger.info(f"Loaded {len(rows)} sample entries from: {path}")
    return rows

def _build_manifest(
    samples: Iterable[Tuple[str, str, str]],
    assemblies_dir: Path,
    logger: logging.Logger,
) -> List[SampleRecord]:
    manifest: List[SampleRecord] = []
    for sample, isolate, specimen in samples:
        genome = assemblies_dir / f"{sample}.fasta"
        exists = genome.exists()
        notes = "" if exists else "FASTA not found"
        if not exists:
            logger.warning(f"Missing genome for sample '{sample}': {genome}")
        manifest.append(SampleRecord(sample, isolate, specimen, genome, exists, notes))
    return manifest

def _write_manifest(manifest: List[SampleRecord], out_path: Path, logger: logging.Logger) -> None:
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "isolate", "specimen", "genome_path", "exists", "notes"])
        for rec in manifest:
            writer.writerow([
                rec.sample,
                rec.isolate,
                rec.specimen,
                str(rec.genome_path),
                "yes" if rec.exists else "no",
                rec.notes,
            ])
    logger.debug(f"Wrote manifest: {out_path}")

def _write_run_config(cfg: Dict, out_path: Path, logger: logging.Logger) -> None:
    out_path.write_text(json.dumps(cfg, indent=2))
    logger.debug(f"Saved run config: {out_path}")

# --- FASTA parsing (pure Python) ---
def _iter_fasta_records(fasta: Path):
    name = None
    seq_chunks: List[str] = []
    with fasta.open("r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, "".join(seq_chunks)

def _first_fasta_id(fa: Path) -> str:
    with fa.open("r") as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    return fa.stem  # fallback

def _prepend_reference_to_manifest(
    manifest: List[SampleRecord],
    reference_fa: Path,
    sample_id: Optional[str],
    logger: logging.Logger,
) -> List[SampleRecord]:
    """Ensure the reference is row 1; if it already exists (by ID or path), move it to front."""
    sid = sample_id or _first_fasta_id(reference_fa)
    ref_path = reference_fa.resolve()
    # find existing by sample ID or by absolute path
    idx = None
    for i, r in enumerate(manifest):
        try:
            if r.sample == sid or r.genome_path.resolve() == ref_path:
                idx = i
                break
        except Exception:
            if r.sample == sid:
                idx = i
                break

    if idx is not None:
        # move existing record to front; don't duplicate
        rec = manifest.pop(idx)
        logger.info(f"Reference '{sid}' already present; moved to first row.")
        return [rec] + manifest

    # otherwise create a new record and prepend
    if not reference_fa.exists():
        logger.error(f"Reference FASTA not found: {reference_fa}")
        return manifest

    ref_rec = SampleRecord(
        sample=sid,
        isolate="reference",
        specimen="reference",
        genome_path=reference_fa,
        exists=True,
        notes=""
    )
    logger.info(f"Reference sample ID: {sid}")
    return [ref_rec] + manifest

# --- MLST (senterica) --- (senterica) ---
def run_mlst(
    genome: Path,
    outdir: Path,
    mlst_db_dir: Path,
    tmpdir: Path,
    logger: logging.Logger,
    scheme: str,
    env_name: Optional[str] = None,
) -> str:
    _ensure_dir(outdir)
    cmd = [
        "mlst.py",
        "-i", str(genome),
        "-o", str(outdir),
        "-s", scheme,               # <- use CAPS scheme
        "-x",
        "-p", str(mlst_db_dir),
        "-t", str(tmpdir),
        "-q",
    ]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (outdir / "mlst.stdout.txt").write_text(out)
    if rc != 0:
        logger.warning(f"MLST failed for {genome.name}")
        return "NA"
    results = outdir / "results.txt"
    if not results.exists():
        logger.warning(f"MLST results.txt missing for {genome.name}")
        return "NA"
    try:
        for line in results.read_text().splitlines():
            if "Sequence Type" in line or line.lower().startswith("st\t"):
                parts = line.replace(":", "\t").split("\t")
                for i, tok in enumerate(parts):
                    if tok.strip().lower() in {"st", "sequence type"} and i + 1 < len(parts):
                        val = parts[i + 1].strip()
                        try:
                            num = int(val)
                            val = "ST" + str(num)
                        except ValueError:
                            # non-numeric -> leave as-is
                            pass
                        return val if val else "NA"
        return "NA"
    except Exception as e:
        logger.warning(f"Failed to parse MLST results for {genome.name}: {e}")
        return "NA"

# --- cgMLST (Salmonella scheme) ---
def run_cgmlst(
    genome: Path,
    outdir: Path,
    cgmlst_db_dir: Path,
    tmpdir: Path,
    logger: logging.Logger,
    scheme: str,
    env_name: Optional[str] = None,
) -> str:
    _ensure_dir(outdir)
    cmd = [
        "cgMLST.py",
        "-i", str(genome),
        "-o", str(outdir),
        "-s", scheme,               # <- use CAPS scheme
        "-db", str(cgmlst_db_dir),
        "-t", str(tmpdir),
        "-q",
    ]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (outdir / "cgmlst.stdout.txt").write_text(out)
    if rc != 0:
        logger.warning(f"cgMLST failed for {genome.name}")
        return "NA"

    summ = outdir / f"{scheme}_summary.txt"
    if not summ.exists():
        logger.warning(f"cgMLST summary not found for {genome.name}")
        return "NA"

    try:
        rows = [r.split("\t") for r in summ.read_text().splitlines() if r.strip()]
        if not rows:
            return "NA"
        header = [h.strip().lower() for h in rows[0]]
        data = rows[1] if len(rows) > 1 else None
        if not data:
            return "NA"
        for key in ["cgst", "cgmlst_st", "cgmlst st", "cgst id", "cg-mlst st"]:
            if key in header:
                idx = header.index(key)
                return data[idx].strip() or "NA"
        if len(data) >= 5:
            return data[4].strip() or "NA"
        return "NA"
    except Exception as e:
        logger.warning(f"Failed to parse cgMLST summary for {genome.name}: {e}")
        return "NA"

def _parse_pairs_header(row: List[str]) -> Tuple[int, int, Optional[int]]:
    if not row:
        return 0, 1, 2
    lower = [c.strip().lower() for c in row]
    f1 = {"ref", "query1", "genome1", "file1", "name1", "a", "input1"}
    f2 = {"query", "query2", "genome2", "file2", "name2", "b", "input2"}
    ani = {"ani", "ani_cluster", "ani_estimate"}

    def find(keys): 
        for k in keys:
            if k in lower: return lower.index(k)
        return None
    
    i1 = find(f1); i2 = find(f2); i_ani = find(ani)
    if i1 is None or i2 is None:
        return 0, 1, 2
    return i1, i2, (i_ani if i_ani is not None else 2)

def write_skani_matrix(pairs_tsv: Path, name_to_sample: Dict[str, str], sample_order: List[str], out_matrix_tsv: Path, logger: logging.Logger) -> None:
    idx = {s: i for i, s in enumerate(sample_order)}
    N = len(sample_order)
    mat = [["-" for _ in range(N)] for __ in range(N)]
    if not pairs_tsv.exists():
        with out_matrix_tsv.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample"] + sample_order)
            for s in sample_order:
                w.writerow([s] + ["-" for _ in range(N)])
        return
    rows = [ln.strip().split("\t") for ln in pairs_tsv.read_text().splitlines() if ln.strip()]
    has_header = bool(rows and any(tok.lower() in {"ani", "af", "ref", "query", "genome1", "genome2", "file1", "file2"} for tok in rows[0]))
    i1, i2, i_ani = _parse_pairs_header(rows[0] if has_header else [])
    for r in (rows[1:] if has_header else rows):
        if len(r) <= max(i1, i2, i_ani):
            continue
        a_name = Path(r[i1]).name; b_name = Path(r[i2]).name
        a = name_to_sample.get(a_name) or Path(a_name).stem
        b = name_to_sample.get(b_name) or Path(b_name).stem
        if a not in idx or b not in idx:
            continue
        ani_val = r[i_ani].strip() or "NA"
        ia, ib = idx[a], idx[b]
        mat[ia][ib] = ani_val; mat[ib][ia] = ani_val
    with out_matrix_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + sample_order)
        for s in sample_order:
            w.writerow([s] + mat[idx[s]])

# --- seq stats ---
def write_assembly_summary(manifest: List[SampleRecord], out_summary_tsv: Path, logger: logging.Logger) -> None:
    """
    Compute per-sample assembly stats directly from FASTAs:
      columns: sample, num_seqs, sum_len, min_len, max_len, gc
      - gc is length-weighted across all contigs (total_GC / total_len * 100)
    """
    header = ["sample", "num_seqs", "sum_len", "min_len", "max_len", "gc"]
    with out_summary_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for rec in manifest:
            if not rec.exists:
                w.writerow([rec.sample, 0, 0, 0, 0, "NA"])
                continue
            n = 0
            total_len = 0
            min_len = None
            max_len = 0
            gc_bases = 0
            try:
                for _name, seq in _iter_fasta_records(rec.genome_path):
                    L = len(seq)
                    n += 1
                    total_len += L
                    min_len = L if (min_len is None or L < min_len) else min_len
                    if L > max_len:
                        max_len = L
                    s = seq.upper()
                    gc_bases += s.count("G") + s.count("C")
                gc_pct = f"{(100.0 * gc_bases / total_len):.2f}" if total_len > 0 else "NA"
                w.writerow([rec.sample, n, total_len, (min_len or 0), max_len, gc_pct])
            except Exception as e:
                logger.warning(f"Failed reading {rec.genome_path}: {e}")
                w.writerow([rec.sample, 0, 0, 0, 0, "NA"])
    logger.info(f"Wrote assembly summary: {out_summary_tsv}")

def _parse_tsv_rows(path: Path) -> List[List[str]]:
    return [ln.strip().split("\t") for ln in path.read_text().splitlines() if ln.strip()]

def run_skani_pair(reference: Path, query: Path, out_tsv: Path, logger: logging.Logger, env_name: Optional[str] = None) -> Tuple[str, str]:
    """
    Return (ANI, AF) as strings. 'NA' on failure.
    We call: skani dist -r REF -q QUERY -o out_tsv
    """
    cmd = ["skani", "dist", "-r", str(reference), "-q", str(query), "-o", str(out_tsv)]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (out_tsv.parent / f"{query.stem}.skani.stdout.txt").write_text(out)
    if rc != 0 or not out_tsv.exists():
        return "NA", "NA"
    try:
        rows = _parse_tsv_rows(out_tsv)
        if not rows:
            return "NA", "NA"

        header = [h.strip() for h in rows[0]]
        data = rows[1] if len(rows) > 1 else rows[0]
        hmap = {h.lower(): i for i, h in enumerate(header)}

        def pick(idx):
            return (data[idx].strip() if idx is not None and idx < len(data) and data[idx].strip() != "" else None)

        # ANI candidates
        ani = None
        for key in ("ani", "average_nucleotide_identity"):
            if key in hmap:
                ani = pick(hmap[key]); break
        if ani is None and len(data) >= 3:
            ani = data[2].strip()

        # AF candidates
        af = None
        if "af" in hmap:
            af = pick(hmap["af"])
        else:
            q = pick(hmap.get("align_fraction_query"))
            r = pick(hmap.get("align_fraction_ref"))
            # choose conservative minimum if both present; else whichever exists
            if q and r:
                try:
                    af = str(min(float(q), float(r)))
                except Exception:
                    af = q
            else:
                af = q or r

        return (ani or "NA"), (af or "NA")
    except Exception:
        return "NA", "NA"

def run_fastani_pair(reference: Path, query: Path, out_tsv: Path, logger: logging.Logger, env_name: Optional[str] = None) -> str:
    """
    Return ANI as string. 'NA' on failure.
    fastANI -r REF -q QUERY -o out_tsv
    Output columns typically: query, ref, ANI, ...
    """
    cmd = ["fastANI", "-r", str(reference), "-q", str(query), "-o", str(out_tsv)]
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    (out_tsv.parent / f"{query.stem}.fastani.stdout.txt").write_text(out)
    if rc != 0 or not out_tsv.exists():
        return "NA"
    try:
        rows = _parse_tsv_rows(out_tsv)
        if not rows:
            return "NA"
        # Usually a single row; third column is ANI
        row = rows[0]
        return row[2] if len(row) >= 3 else "NA"
    except Exception:
        return "NA"

def write_reference_ani_table(manifest: List[SampleRecord],
                              reference_index: int,
                              ani_tool: str,
                              out_tsv: Path,
                              workdir: Path,
                              logger: logging.Logger,
                              skani_env: Optional[str] = None,
                              fastani_env: Optional[str] = None) -> None:
    """
    For each sample vs manifest[reference_index], compute ANI.
    Reference row gets a dash '-'.
    Output columns: sample, ANI, AF, tool
    """
    ref = manifest[reference_index]
    workdir.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "ANI", "AF", "tool"])
        # reference row
        w.writerow([ref.sample, "-", "-", ani_tool if ani_tool != "auto" else "auto"])
        # others
        for i, rec in enumerate(manifest):
            if i == reference_index:
                continue
            if not rec.exists:
                w.writerow([rec.sample, "NA", "NA", ani_tool])
                continue
            # temp per-pair outfile
            pair_out = workdir / f"{ref.sample}__vs__{rec.sample}.tsv"
            ani = af = "NA"
            tool_used = "NA"
            if ani_tool in ("auto", "skani"):
                ani, af = run_skani_pair(ref.genome_path, rec.genome_path, pair_out, logger, env_name=skani_env)
                tool_used = "skani"
                if ani_tool == "auto" and ani == "NA":
                    # fallback to fastANI
                    ani = run_fastani_pair(ref.genome_path, rec.genome_path, pair_out, logger, env_name=fastani_env)
                    af = "NA"
                    tool_used = "fastANI" if ani != "NA" else "auto"
            else:
                ani = run_fastani_pair(ref.genome_path, rec.genome_path, pair_out, logger, env_name=fastani_env)
                af = "NA"
                tool_used = "fastANI"
            w.writerow([rec.sample, ani, af, tool_used])
    logger.info(f"Wrote ANI vs reference: {out_tsv}")

def write_skani_pairs_with_samples(
    pairs_tsv: Path,
    name_to_sample: Dict[str, str],
    out_tsv: Path,
    logger: logging.Logger,
) -> None:
    """Rewrite skani pairs file to sample-labeled table: sample1, sample2, ani, af."""
    if not pairs_tsv.exists():
        logger.warning(f"skani pairs not found: {pairs_tsv}")
        out_tsv.write_text("sample1\tsample2\tani\taf\n")
        return
    rows = [ln.strip().split("\t") for ln in pairs_tsv.read_text().splitlines() if ln.strip()]
    if not rows:
        out_tsv.write_text("sample1\tsample2\tani\taf\n")
        return

    # header detection (reuse _parse_pairs_header)
    has_header = any(tok.lower() in {"ani","af","ref","query","genome1","genome2","file1","file2"} for tok in rows[0])
    i1, i2, i_ani = _parse_pairs_header(rows[0] if has_header else [])

    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample1", "sample2", "ani", "af"])
        it = rows[1:] if has_header else rows
        for r in it:
            if len(r) <= max(i1, i2, i_ani):
                continue
            a_name = Path(r[i1]).name
            b_name = Path(r[i2]).name
            a = name_to_sample.get(a_name) or Path(a_name).stem
            b = name_to_sample.get(b_name) or Path(b_name).stem
            ani = r[i_ani].strip() if len(r) > i_ani else "NA"
            # AF may be absent in some formats
            af = "NA"
            if "af" in [c.lower() for c in rows[0]] if has_header else False:
                j_af = [c.lower() for c in rows[0]].index("af")
                if len(r) > j_af:
                    af = r[j_af].strip() or "NA"
            elif len(r) >= 4:
                af = r[3].strip() or "NA"
            w.writerow([a, b, ani or "NA", af or "NA"])

# ----------------------------
# Orchestrator per sample
# ----------------------------
def _make_paths(accession: str, outdir: Path, tmpdir: Path) -> Dict[str, Path]:
    p = {
        "outdir": outdir,
        "serotype_root": outdir / "serotype",
        "mlst_root": outdir / "mlst",
        "cgmlst_root": outdir / "cgMLST",
        "tmpdir": tmpdir,

        "gc_tsv": outdir / f"{accession}_GC.tsv",
        "serotypes_tsv": outdir / f"{accession}_serotypes.tsv",
        "mlst_tsv": outdir / f"{accession}_mlst.tsv",
        "cgmlst_tsv": outdir / f"{accession}_cgmlst.tsv",
        "assembly_summary_tsv": outdir / f"{accession}_assembly_summary.tsv",

        "ani_tsv": outdir / f"{accession}_ani.tsv",
        "ani_pairs_dir": outdir / "ani",

        # skani
        "skani_root": outdir / "ani",
        "skani_pairs": outdir / "ani" / f"{accession}_skani_pairs.tsv",
        "skani_pairs_samples": outdir / "ani" / f"{accession}_skani_pairs_samples.tsv",
        "skani_matrix": outdir / "ani" / f"{accession}_skani_matrix.tsv",

        # SKA2
        "ska_root": outdir / "ska",
        "ska_filelist": outdir / "ska" / f"{accession}_ska.filelist.txt",
        "ska_skf": outdir / "ska" / f"{accession}.skf",
        "ska_distances": outdir / "ska" / f"{accession}_ska.distances.tsv",
        "ska_vs_ref": outdir / "ska" / f"{accession}_ska_vs_reference.tsv",
        "ska_matrix_snps": outdir / "ska" / f"{accession}_ska_matrix.snps.tsv",
        "ska_matrix_mash": outdir / "ska" / f"{accession}_ska_matrix.mash.tsv",
    }
    for k in ("serotype_root", "mlst_root", "cgmlst_root", "ani_pairs_dir", "ska_root"):
        p[k].mkdir(parents=True, exist_ok=True)
    return p

def _load_prev_map_2col(path: Path) -> Dict[str, Tuple[str, str]]:
    # For serotypes.tsv: sample, serotype, formula
    if not path.exists(): return {}
    d = {}
    with path.open() as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r, None)
        for row in r:
            if len(row) >= 3:
                d[row[0]] = (row[1], row[2])
    return d

def _load_prev_map_1col(path: Path) -> Dict[str, str]:
    # For mlst.tsv (sample, ST) and cgmlst.tsv (sample, cgST)
    if not path.exists(): return {}
    d = {}
    with path.open() as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r, None)
        for row in r:
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

def _print_aligned_table(headers: List[str], rows: List[List[str]], pad: int = 2) -> None:
    """Print a left-aligned, space-padded table sized to the widest cell per column."""
    # ensure strings
    rows = [[str(c) for c in r] for r in rows]
    headers = [str(h) for h in headers]
    ncols = len(headers)
    # compute column widths
    widths = [len(headers[i]) for i in range(ncols)]
    for r in rows:
        for i in range(ncols):
            if i < len(r):
                w = len(r[i])
                if w > widths[i]:
                    widths[i] = w
    # build format string
    sep = " " * pad
    fmt = sep.join(f"{{:<{w}}}" for w in widths)

    # print header + underline + rows
    print(fmt.format(*headers))
    print(fmt.format(*("-" * w for w in widths)))
    for r in rows:
        # pad missing cells if any
        r += [""] * (ncols - len(r))
        print(fmt.format(*r))

def prepare_ska_filelist_preflight(
    manifest: List[SampleRecord],
    ska_root: Path,
    out_prefix: str,
    logger: logging.Logger,
) -> Tuple[Path, int]:
    """
    Write <ska_root>/<out_prefix>.filelist.txt with absolute genome paths for all
    existing samples from manifest. Returns (file_list_path, n_entries).
    """
    ska_root.mkdir(parents=True, exist_ok=True)
    out_name = Path(out_prefix).name
    file_list = ska_root / f"{out_name}.filelist.txt"

    n = 0
    with file_list.open("w") as fh:
        for rec in manifest:
            if not rec.exists:
                continue
            fh.write(str(rec.genome_path.resolve()) + "\n")
            n += 1

    logger.info(f"SKA preflight: wrote {n} genomes -> {file_list}")
    return file_list, n

def _load_manifest_tsv(path: Path, logger: logging.Logger) -> List[SampleRecord]:
    """Read <run>_manifest.tsv back into SampleRecord list."""
    out: List[SampleRecord] = []
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open() as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r, None) or []
        h = {k.strip().lower(): i for i, k in enumerate(header)}
        for row in r:
            if not row: 
                continue
            sample   = row[h.get("sample", 0)]
            isolate  = row[h.get("isolate", 1)] if "isolate" in h else ""
            specimen = row[h.get("specimen", 2)] if "specimen" in h else ""
            date     = row[h.get("date", 3)] if "date" in h else ""
            gp_idx   = h.get("genome_path")
            gpath    = Path(row[gp_idx]) if gp_idx is not None and gp_idx < len(row) else Path("")
            ex_idx   = h.get("exists")
            exists   = (row[ex_idx].strip().lower() in {"yes", "true", "1"}) if ex_idx is not None else gpath.exists()
            notes    = row[h.get("notes", len(row)-1)] if "notes" in h else ""
            out.append(SampleRecord(sample, isolate, specimen, gpath, exists, notes))
    logger.info(f"Loaded existing manifest for resume: {path} ({len(out)} rows)")
    return out

# ----------------------------
# CLI
# ----------------------------
def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="type-genomes",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent(
            """
            Set up a typing run from assembled genomes and run core typing.

            Inputs: sample sheet (sample, isolate, specimen, date) and <sample>.fasta in assemblies dir.
            Outputs: manifest, run config, tool subdirs, and consolidated TSVs in outdir.
            """
        ),
    )
    req = p.add_argument_group("Required arguments")
    opt = p.add_argument_group("Optional arguments")

    req.add_argument("-n", "--name", dest="tech_name",
                help="Analyst / technician name shown in the report header.")

    req.add_argument("--run-name", help="Run name identifier (e.g., 250701_13_SEQ10_LMO)")

    req.add_argument("-a", "--accession", dest="accession",
                help="Accession / Case ID / output reports prefix.")

    req.add_argument("--organism", help="Organism label name). Check spelling with --list-organisms.")
    req.add_argument("--assemblies-dir", type=Path, 
    help="<sample>.fasta directory.")
    req.add_argument("--sample-sheet", type=Path, help="4-column TSV. Default: try common Dropbox locations.")
    opt.add_argument("--outdir", type=Path, help="Output directory. Default: ./bactipipe_relate_output/Typing_<run-name>/<accession>")
    opt.add_argument("--tmpdir", type=Path, help="Temporary work directory. Default: <outdir>/tmp")
    opt.add_argument("--reference", type=Path, help="Optional reference genome (treated as first sample).")
    opt.add_argument("--reference-id", help="Optional sample ID for the reference. If omitted, uses the first FASTA header ID.")
    opt.add_argument("--ani-tool", choices=["auto", "skani", "fastani"], default="auto",
                   help="ANI vs reference tool (auto prefers skani, else fastANI).")
    opt.add_argument("--db-dir", type=Path, help="DB root directory. Defaults to $BACTIPIPE_DB_DIR")
    opt.add_argument("--threads", type=int, default=4, help="Threads for serotyping (SeqSero2).")
    opt.add_argument("--list-organisms", action="store_true", help="List supported organisms/schemes and exit")
    opt.add_argument("--no-ska", action="store_true", help="Disable SKA split-kmer relatedness (enabled by default).")
    opt.add_argument("-v", "--verbose", action="store_true", help="More console detail (DEBUG).")
    opt.add_argument(
    "--resume",
    action="store_true",
    help="Skip heavy steps (serotype/MLST/cgMLST/ANI/SKA) when final outputs already exist; parse and reuse them.",
    )
    opt.add_argument("--title", default="Strain Relatedness Analysis",
                help="Report title (used in PDF header).")

    # Hidden env flags (exception tools only). Defaults can come from env vars; users can still pass them explicitly.
    opt.add_argument("--seqsero2-env",
                   default=os.environ.get("BACTIPIPE_ENV_SEQSERO2") or "genepid",
                   help=argparse.SUPPRESS)
    opt.add_argument("--kleborate-env",
                   default=os.environ.get("BACTIPIPE_ENV_KLEBORATE") or "bactipipe",
                   help=argparse.SUPPRESS)    
    opt.add_argument("--cge-env",
                   default=os.environ.get("BACTIPIPE_ENV_CGE") or "genepid",
                   help=argparse.SUPPRESS)
<<<<<<< HEAD
=======
    opt.add_argument("--cge-env",
                   default=os.environ.get("BACTIPIPE_ENV_CGE") or "genepid",
                   help=argparse.SUPPRESS)
>>>>>>> 1ea388f (removed false AIDA (AIDA-like) records)

    return p

def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_argparser()
    args = ap.parse_args(argv)

    # --- short-circuit: just list organisms and exit ---
    if args.list_organisms:
        headers = ["organism", "uses_mlst", "mlst_scheme", "uses_cgmlst", "cgmlst_scheme", "serotyper"]
        rows = supported_rows()  # you already normalized Nones inside this
        _print_aligned_table(headers, rows, pad=2)
        return 0

    # --- validate required args only for real runs ---
    missing = []
    if not args.tech_name:      missing.append("--name")
    if not args.assemblies_dir:  missing.append("--assemblies-dir")
    if not args.sample_sheet: missing.append("--sample-sheet")
    if not args.organism:    missing.append("--organism")
    if not args.run_name:     missing.append("--run-name")
    if not args.accession:    missing.append("--accession")

    if missing:
        sys.stderr.write("ERROR: missing required arguments: " + ", ".join(missing) + "\n")
        sys.exit(2)


    # Resolve defaults
    run_name: str = args.run_name
    accession: str = args.accession
    assemblies_dir = args.assemblies_dir or Path(f"./{run_name}/assemblies/genomes")
    outdir = args.outdir or Path(f"./Typing_{run_name}/{accession}")
    tmpdir = args.tmpdir or (outdir / "tmp")

    # Sample sheet default
    sample_sheet = args.sample_sheet or _default_sample_sheet(run_name)
    if sample_sheet is None:
        ap.error("--sample-sheet not provided and no default found. Set --sample-sheet or BACTIPIPE_SAMPLE_SHEETS_DIR.")

    # DB dir default
    db_dir = args.db_dir
    if db_dir is None:
        db_env = os.environ.get("BACTIPIPE_DB_DIR")
        if db_env:
            db_dir = Path(db_env)
        else:
            ap.error("--db-dir not provided and BACTIPIPE_DB_DIR not set.")


    # Hard-coded subdirs (no CLI overrides)
    serotypefinder_db_dir = db_dir / "serotypefinder_db"
    mlst_db_dir          = db_dir / "mlst_db"
    cgmlst_db_dir        = db_dir / "cgmlstfinder_db"

    org = args.organism.lower()
    caps = CAPS.get(org, {})
    need_mlst   = bool(caps.get("mlst"))
    need_cgmlst = bool(caps.get("cgmlst"))
    mlst_scheme   = caps.get("mlst_scheme") if need_mlst else None
    cgmlst_scheme = caps.get("cgmlst_scheme") if need_cgmlst else None

    # Determine which serotyper this organism uses
    serotyper_label = (caps.get("serotyper") or "").lower()
    # Default (unknown/empty serotyper) falls back to SerotypeFinder for non-Salmonella
    uses_serofinder = (serotyper_label == "serotypefinder") or (
        serotyper_label == "" and org not in {"salmonella"}
    )

    run_cfg = {
        "run_name": run_name,
        "organism": args.organism,
        "assemblies_dir": str(assemblies_dir.resolve()),
        "sample_sheet": str(sample_sheet.resolve()),
        "outdir": str(outdir.resolve()),
        "tmpdir": str(tmpdir.resolve()),
        "reference": str(args.reference.resolve()) if args.reference else None,
        "db_dir": str(db_dir.resolve()),
        "mlst_db_dir": str(mlst_db_dir) if need_mlst else None,
        "cgmlst_db_dir": str(cgmlst_db_dir) if need_cgmlst else None,
        "serotypefinder_db_dir": str(serotypefinder_db_dir) if uses_serofinder else None,
        "threads": args.threads,
        "phase": 2,
        "envs": {
            "seqsero2": args.seqsero2_env,
            "kleborate": args.kleborate_env,
            "mlst": "genepid",
            "cgmlst": "genepid",
        },       
    }

    ani_env_skani   = getattr(args, "skani_env", None)    # usually None
    ani_env_fastani = getattr(args, "fastani_env", None)  # usually None

    # Init logging
    logger = _init_logging(outdir, verbose=args.verbose)

    config_path = outdir / f"{accession}_run_config.json"
    _write_run_config(run_cfg, config_path, logger)

    # Preflight: check required tools
    logger.info("Preflight: checking required tools")
    ok = True
    def _must(name, env=None):
        ok_local = _check_tool(name, logger, env_name=env)
        logger.info(f"  - {name:<18} {'ok' if ok_local else 'MISSING'}" + (f"  [env={env}]" if env else ""))
        return ok_local

    # serotyper (by organism)
    serotyper_label = (caps.get("serotyper") or "").lower()
    if args.organism.lower() == "salmonella":
        ok &= _must("SeqSero2_package.py", env=run_cfg["envs"].get("seqsero2"))
    elif serotyper_label == "kleborate":
        ok &= _must("kleborate", env=run_cfg["envs"].get("kleborate"))
    else:
        ok &= _must("serotypefinder.py", env=run_cfg["envs"].get("serotypefinder"))

    # MLST / cgMLST (if used)
    if need_mlst:   ok &= _must("mlst.py",   env=run_cfg["envs"].get("mlst"))
    if need_cgmlst: ok &= _must("cgMLST.py", env=run_cfg["envs"].get("cgmlst"))

    # ANI tools
    has_ref = bool(args.reference)
    ani_mode = args.ani_tool
    if has_ref:
        if args.ani_tool == "skani":
            ok &= _must("skani", env=ani_env_skani)
        elif args.ani_tool == "fastani":
            ok &= _must("fastANI", env=ani_env_fastani)
        else:  # auto
            if not _must("skani", env=ani_env_skani):
                ok &= _must("fastANI", env=ani_env_fastani)
    else:
        ok &= _must("skani", env=ani_env_skani)

    # SKA2 (unless disabled)
    if not args.no_ska:
        ok &= _must("ska")  # ska2 binary

    if not ok:
        logger.error("Preflight failed: missing required tool(s). Exiting.")
        return 2
    logger.info("Preflight ✓\n")

    # layout
    logger.info(f"Run name: {run_name}")
    logger.info(f"Accession: {accession}")
    logger.info(f"Organism: {args.organism}")
    logger.info(f"Assemblies dir: {assemblies_dir}")
    logger.info(f"Sample sheet: {sample_sheet}")
    logger.info(f"Outdir: {outdir}")
    logger.info(f"Tmpdir: {tmpdir}")
    logger.info(f"DB dir: {db_dir}")
    if need_mlst:   logger.info(f"MLST DB: {mlst_db_dir}")
    if need_cgmlst: logger.info(f"cgMLST DB: {cgmlst_db_dir}")
    if uses_serofinder: logger.info(f"SerotypeFinder DB: {serotypefinder_db_dir}")

    # Validate DB paths early (only those required by organism features)
    err = False
    if need_mlst and not mlst_db_dir.exists():
        logger.error(f"MLST DB not found: {mlst_db_dir}"); err = True
    if need_cgmlst and not cgmlst_db_dir.exists():
        logger.error(f"cgMLST DB not found: {cgmlst_db_dir}"); err = True
    if uses_serofinder and not serotypefinder_db_dir.exists():
        logger.error(f"SerotypeFinder DB not found: {serotypefinder_db_dir}"); err = True

    # SerotypeFinder DB is only needed for organisms using it via dispatch
    uses_serofinder = caps.get("serotyper", "").lower() in {"serotypefinder", "sero", "finder"}
    if uses_serofinder and not serotypefinder_db_dir.exists():
        logger.error(f"SerotypeFinder DB not found: {serotypefinder_db_dir}"); err = True
    if err:
        return 2

    if args.reference:
        logger.info(f"Reference: {args.reference} (treated as first sample)")


    _ensure_layout(outdir, create_tmp=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    # Parse samples and build manifest (Phase 1)
    try:
        rows = _read_sample_sheet(sample_sheet, logger)
    except Exception as e:
        logger.exception(f"Failed to read sample sheet: {e}")
        return 2

    manifest = _build_manifest(rows, assemblies_dir, logger)
    if args.reference:
        manifest = _prepend_reference_to_manifest(manifest, args.reference, args.reference_id, logger)
    # If a reference is provided, prepend it as row 1 and run it like any sample
    if args.reference:
        manifest = _prepend_reference_to_manifest(manifest, args.reference, args.reference_id, logger)

    # Persist manifest + run config
    manifest_path = outdir / f"{accession}_manifest.tsv"

    if args.resume and manifest_path.exists():
        manifest = _load_manifest_tsv(manifest_path, logger)
    else:
        # fresh build from sample sheet
        try:
            rows = _read_sample_sheet(sample_sheet, logger)
        except Exception as e:
            logger.exception(f"Failed to read sample sheet: {e}")
            return 2
        manifest = _build_manifest(rows, assemblies_dir, logger)
        if args.reference:
            manifest = _prepend_reference_to_manifest(manifest, args.reference, args.reference_id, logger)
        # persist only on fresh build (resume keeps prior)
        _write_manifest(manifest, manifest_path, logger)

    # --- SKA preflight (fail fast if enabled but unusable) ---
    if not args.no_ska:
        ska_root = outdir / "ska"  # local var (don't use paths yet)
        with _step(logger, "SKA preflight (file list)"):
            ska_filelist, ska_n = prepare_ska_filelist_preflight(
                manifest=manifest,
                ska_root=ska_root,
                out_prefix=f"{accession}_ska",
                logger=logger,
            )
            if ska_n < 2:
                logger.error("SKA requires at least 2 existing genomes. "
                            "Either add more genomes or rerun with --no-ska to skip SKA.")
                return 2
    else:
        ska_root = outdir / "ska"  
        ska_filelist = None

    # Prepare paths for outputs
    paths = _make_paths(accession, outdir, tmpdir)
    if ska_filelist is not None:
        paths["ska_filelist"] = ska_filelist

    # --- Per-sample typing for ALL organisms ---
    prev_sero = _load_prev_map_2col(paths["serotypes_tsv"]) if args.resume else {}
    prev_mlst = _load_prev_map_1col(paths["mlst_tsv"])      if args.resume else {}
    prev_cgst = _load_prev_map_1col(paths["cgmlst_tsv"])    if args.resume else {}

    # (Re)initialize the top-level TSVs
    with paths["serotypes_tsv"].open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow(["sample", "serotype", "formula"])
    with paths["mlst_tsv"].open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow(["sample", "ST"])
    with paths["cgmlst_tsv"].open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow(["sample", "cgST"])


    total = len(manifest)
    for idx, rec in enumerate(manifest, start=1):
        logger.info(f"[{idx}/{total}] Sample: {rec.sample}")
        # Serotyping
        if rec.sample in prev_sero:
            serotype, formula = prev_sero[rec.sample]
        else:
            sero, formula, species = serotype_dispatch(
                organism=args.organism,
                genome=rec.genome_path,
                outdir=paths["serotype_root"] / rec.sample,
                db_root=Path(run_cfg["db_dir"]),
                logger=logger,
                threads=args.threads,
                envs=run_cfg["envs"],
            )
            serotype = sero

        with paths["serotypes_tsv"].open("a", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerow([rec.sample, serotype, formula])

        # MLST (only if enabled for this organism)
        if need_mlst:
            if rec.sample in prev_mlst:
                st = prev_mlst[rec.sample]
            else:
                st = run_mlst(
                    rec.genome_path,
                    paths["mlst_root"] / rec.sample,
                    Path(run_cfg["mlst_db_dir"]),
                    paths["tmpdir"],
                    logger,
                    scheme=CAPS.get(args.organism.lower(), {}).get("mlst_scheme", "senterica"),
                    env_name=run_cfg["envs"].get("mlst"),
                )
            with paths["mlst_tsv"].open("a", newline="") as fh:
                csv.writer(fh, delimiter="\t").writerow([rec.sample, st])

        # cgMLST (only if enabled)
        if need_cgmlst:
            if rec.sample in prev_cgst:
                cgst = prev_cgst[rec.sample]
            else:
                cgst = run_cgmlst(
                    rec.genome_path,
                    paths["cgmlst_root"] / rec.sample,
                    Path(run_cfg["cgmlst_db_dir"]),
                    paths["tmpdir"],
                    logger,
                    scheme=CAPS.get(args.organism.lower(), {}).get("cgmlst_scheme"),
                    env_name=run_cfg["envs"].get("cgmlst"),
                )
            with paths["cgmlst_tsv"].open("a", newline="") as fh:
                csv.writer(fh, delimiter="\t").writerow([rec.sample, cgst])
    # =========================
    # ANI logic (verbose steps)
    # =========================
    if args.reference:
        with _step(logger, "ANI vs reference (skani/fastANI)"):
            write_reference_ani_table(
                manifest=manifest,
                reference_index=0,
                ani_tool=args.ani_tool,
                out_tsv=paths["ani_tsv"],
                workdir=paths["ani_pairs_dir"],
                logger=logger,
                skani_env=ani_env_skani,
                fastani_env=ani_env_fastani,
            )
            logger.info(f"  ANI table: {paths['ani_tsv']}")
    else:
        with _step(logger, "ANI triangle (skani all-vs-all)"):
            # 1) Build (sample_name, fasta_path) list
            samples = [(rec.sample, rec.genome_path) for rec in manifest if rec.exists]
            logger.info(f"  Genomes for triangle: {len(samples)}")
            if len(samples) < 2:
                logger.info("  Not enough genomes (need ≥ 2); skipping triangle.")
            else:
                # 2) Ensure skani_root path exists in paths dict
                paths["skani_root"] = paths.get("skani_root", paths["outdir"] / "skani")

                # 3) Run triangle and write *.pairs / *.pairs_samples
                skani_out = U.run_skani_triangle_and_write(
                    samples=samples,
                    skani_root=paths["skani_root"],
                    out_prefix=f"{accession}_skani",
                    logger=logger,
                    env_name=(args.skani_env if hasattr(args, "skani_env") else None),
                )

                if not skani_out:
                    logger.error("  skani triangle failed; no output produced.")
                else:
                    # 4) Publish outputs into paths for the final writer
                    if "pairs" in skani_out:
                        paths["skani_pairs"] = skani_out["pairs"]
                        logger.info(f"  skani pairs (files):   {paths['skani_pairs']}")
                    if "pairs_samples" in skani_out:
                        paths["skani_pairs_samples"] = skani_out["pairs_samples"]
                        logger.info(f"  skani pairs (samples): {paths['skani_pairs_samples']}")
                    if "triangle" in skani_out:
                        logger.info(f"  skani triangle (raw):  {skani_out['triangle']}")
    # === SKA2 split-kmer relatedness (always on unless disabled) ===
    if not args.no_ska:
        logger.info("SKA2 build + distance ...")
        genomes = [rec.genome_path for rec in manifest if rec.exists]
        sample_records = [(rec.sample, rec.genome_path) for rec in manifest if rec.exists] 
        if len(genomes) < 2:
            logger.info("SKA: fewer than 2 genomes present; skipping.")
        else:
            dist_path, ska_stdout = ska2.run_ska_allpairs(
                sample_records=sample_records,
                ska_root=paths["ska_root"],
                out_prefix=f"{accession}_ska",
                k=31,
                threads=None,
            )
            (paths["ska_root"] / "ska.stdout.txt").write_text(ska_stdout)

            if dist_path and dist_path.exists():
                # Optional vs-reference table
                if args.reference:
                    ska2.write_ska_vs_reference(
                        distances_tsv=dist_path,
                        reference_sample=manifest[0].sample,   # reference is first row
                        out_tsv=paths["ska_vs_ref"],
                        logger=logger,
                        sample_order=[rec.sample for rec in manifest if rec.exists],
                    )
                    logger.info(f"  SKA vs-ref: {paths['ska_vs_ref']}")

                # Matrices
                sample_order = [rec.sample for rec in manifest if rec.exists]
                made_snps = ska2.write_ska_matrix(paths["ska_distances"], sample_order, paths["ska_matrix_snps"], logger, metric="snps")
                made_mm   = ska2.write_ska_matrix(paths["ska_distances"], sample_order, paths["ska_matrix_mash"], logger, metric="mismatch")
                if made_snps:
                    logger.info(f"  SKA SNPs matrix: {paths['ska_matrix_snps']}")
                if made_mm:
                    logger.info(f"  SKA mismatch matrix: {paths['ska_matrix_mash']}")
            else:
                logger.warning("SKA: distance failed or output missing.")
    else:
        logger.info("SKA disabled by --no-ska.")

    # Write assembly summary
    write_assembly_summary(manifest, paths["assembly_summary_tsv"], logger)

    with _step(logger, "Final summary table"):
        ska_ref_available = bool(args.reference) and (not args.no_ska) and paths["ska_vs_ref"].exists()
        final_path = U.write_final_summary(
            accession=accession,
            manifest=manifest,
            paths=paths,
            logger=logger,
            reference_provided=bool(args.reference),
        )

    # Final summary to console
    total = len(manifest)
    missing = sum(1 for r in manifest if not r.exists)
    logger.info("")
    logger.info("====== Pipeline summary ======")
    logger.info(f"Samples: {total}; Genomes missing: {missing}")
    logger.info(f"Serotypes TSV: {paths['serotypes_tsv']}")
    logger.info(f"MLST TSV:     {paths['mlst_tsv']}")
    logger.info(f"cgMLST TSV:   {paths['cgmlst_tsv']}")
    logger.info(f"GC TSV:       {paths['gc_tsv']}")
    logger.info(f"Assembly summary TSV: {paths['assembly_summary_tsv']}")
    if args.reference:
        logger.info(f"ANI vs reference: {paths['ani_tsv']}")
    
    # === PDF report ===
    with _step(logger, "Render PDF report"):
        header_block = U.make_report_header_block(
            title=args.title,
            run_name=args.run_name,
            tech_name=args.tech_name,
            accession=args.accession or "",
        )
        used = U.detect_used_tools(
            paths,
            logger,
            env_overrides={
                "seqsero2": getattr(args, "seqsero2_env", None),
                "cge":      getattr(args, "cge_env", None),
                "kleborate":getattr(args, "kleborate_env", None),
                "skani":    getattr(args, "skani_env", None),
                "ska":      getattr(args, "ska_env", None),
            },
        )
        versions = U.collect_tool_versions(used, logger)

        final_tsv = outdir / f"{accession}_relate.tsv"
        out_pdf   = outdir / f"{accession}_relate.pdf"

        U.render_pdf_type_genomes_full(
            final_tsv=final_tsv,
            out_pdf=out_pdf,
            title=args.title,
            header_text=header_block,
            tool_versions=versions,
            run_name=run_name,
            accession=accession,
        )

    logger.info("Pipeline complete.")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
