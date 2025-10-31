from __future__ import annotations
from typing import Dict, Tuple, List, Optional
import csv
import logging
from pathlib import Path
import subprocess
from bactipipe.scripts.utils import _fmt_int_commas

def _parse_ska_distances(distances_tsv: Path, logger: logging.Logger) -> Dict[Tuple[str,str], Dict[str,str]]:
    """
    Parse SKA2 'distance' TSV into a dict:
      {(s1, s2): {"snps": "<int-like>", "mismatch": "<float-like>"}}
    Recognizes headers: Sample1/Sample2 and Distance/Mismatches (SKA2),
    and legacy synonyms snps/snp_distance, mash/mash_distance.
    Provides symmetric lookups and 0 on the diagonal.
    """
    pairs: Dict[Tuple[str,str], Dict[str,str]] = {}
    if not distances_tsv.exists():
        logger.warning(f"SKA distances file not found: {distances_tsv}")
        return pairs

    with distances_tsv.open() as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if any(c.strip() for c in r)]
    if not rows:
        return pairs

    # Header detection
    header = [h.strip() for h in rows[0]]
    header_lower = [h.lower() for h in header]
    has_header = any(h in header_lower for h in ("sample1","sample2","distance","mismatches","snps","mash","snp_distance","mash_distance"))
    start = 1 if has_header else 0

    hmap = {h.lower(): i for i, h in enumerate(header)} if has_header else {}

    # Column resolvers (prefer SKA2 names, then legacy)
    def col_idx(*names: str) -> Optional[int]:
        for name in names:
            i = hmap.get(name)
            if i is not None:
                return i
        return None

    s1_idx = col_idx("sample1", "sample_a", "name1") if has_header else 0
    s2_idx = col_idx("sample2", "sample_b", "name2") if has_header else 1
    snp_idx = col_idx("distance", "snps", "snp_distance")
    mm_idx  = col_idx("mismatches", "mash", "mash_distance")

    def is_num(s: str) -> bool:
        try:
            float(s)
            return True
        except Exception:
            return False

    for r in rows[start:]:
        if len(r) < 3:
            continue

        if has_header:
            s1 = r[s1_idx].strip() if s1_idx is not None and s1_idx < len(r) else ""
            s2 = r[s2_idx].strip() if s2_idx is not None and s2_idx < len(r) else ""
            snps = (r[snp_idx].strip() if snp_idx is not None and snp_idx < len(r) else "NA")
            mismatch = (r[mm_idx].strip()  if mm_idx  is not None and mm_idx  < len(r) else "NA")
        else:
            # No header fallback: s1, s2, snps[, mismatch]
            s1 = r[0].strip()
            s2 = r[1].strip()
            snps = r[2].strip() if len(r) >= 3 and is_num(r[2]) else "NA"
            mismatch = r[3].strip() if len(r) >= 4 and is_num(r[3]) else "NA"

        if not s1 or not s2:
            continue

        pairs[(s1, s2)] = {"snps": snps, "mismatch": mismatch}
        pairs[(s2, s1)] = {"snps": snps, "mismatch": mismatch}
        pairs.setdefault((s1, s1), {"snps": "0", "mismatch": "0"})
        pairs.setdefault((s2, s2), {"snps": "0", "mismatch": "0"})

    return pairs


def write_ska_vs_reference(
    distances_tsv: Path,
    reference_sample: str,
    out_tsv: Path,
    logger: logging.Logger,
    sample_order: Optional[List[str]] = None,
) -> bool:
    """
    Write a 3-column table of SKA distances vs reference:
      sample, snps, mismatch
    """
    pairs = _parse_ska_distances(distances_tsv, logger)
    if not pairs:
        return False

    if sample_order:
        order = [s for s in sample_order if s != reference_sample]
    else:
        seen = {a for (a, _) in pairs.keys()} | {b for (_, b) in pairs.keys()}
        order = sorted([s for s in seen if s != reference_sample])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "snps", "mismatch"])
        for s in order:
            d = pairs.get((reference_sample, s), {"snps": "NA", "mismatch": "NA"})
            snps = _fmt_int_commas(d.get("snps", "NA"))        # ← add commas
            mm   = d.get("mismatch", "NA")                     # keep as-is (fraction)
            w.writerow([s, snps, mm])
    return True


def write_ska_matrix(
    distances_tsv: Path,
    sample_order: List[str],
    out_tsv: Path,
    logger: logging.Logger,
    metric: str = "snps",   # accepts 'snps' or 'mismatch'
) -> bool:
    metric = metric.lower()
    if metric not in {"snps", "mismatch"}:
        raise ValueError("metric must be 'snps' or 'mismatch'")

    pairs = _parse_ska_distances(distances_tsv, logger)
    if not pairs:
        return False

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + sample_order)
        for a in sample_order:
            row = [a]
            for b in sample_order:
                if a == b:
                    cell = "0"
                else:
                    d = pairs.get((a, b))
                    cell = d.get(metric, "NA") if d else "NA"
                if metric == "snps":
                    cell = _fmt_int_commas(cell)               # ← add commas
                row.append(cell)
            w.writerow(row)
    return True


def _run(cmd: List[str], cwd: Optional[Path] = None) -> Tuple[int, str]:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return p.returncode, p.stdout

def write_filelist(sample_records: List[Tuple[str, Path]], filelist: Path) -> Path:
    """
    Write SKA2 build file list:
      <sample_name>\t</abs/path/to/genome.fa>
    One sample per line.
    """
    filelist.parent.mkdir(parents=True, exist_ok=True)
    with filelist.open("w") as fh:
        for sample, g in sample_records:
            fh.write(f"{sample}\t{g.resolve()}\n")
    return filelist

def ska2_build(prefix: Path, filelist: Path, k: int = 31, single_strand: bool = False) -> Tuple[int, str, Path]:
    """
    Build an SKF: ska build -o <prefix> -k <k> -f <filelist> [--single-strand]
    Returns (rc, stdout, skf_path).
    """
    cmd = ["ska", "build", "-o", str(prefix), "-k", str(k), "-f", str(filelist)]
    if single_strand:
        # '--single-strand' is valid; place it right after 'build' if you prefer
        cmd.insert(2, "--single-strand")
    rc, out = _run(cmd)
    skf = prefix if prefix.suffix == ".skf" else prefix.with_suffix(".skf")
    return rc, out, skf

def ska2_distance(
    skf: Path,
    out_tsv: Path,
    threads: Optional[int] = None,
    min_freq: Optional[float] = None,
    allow_ambiguous: bool = False,
) -> Tuple[int, str]:
    """
    Calculate all-pairs distances from an SKF:
      ska distance <SKF_FILE> -o <out.tsv> [--threads N] [--min-freq F] [--allow-ambiguous]
    Returns (rc, stdout).
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["ska", "distance", str(skf), "-o", str(out_tsv)]
    if threads is not None:
        cmd += ["--threads", str(threads)]
    if min_freq is not None:
        cmd += ["--min-freq", str(min_freq)]
    if allow_ambiguous:
        cmd += ["--allow-ambiguous"]
    return _run(cmd)

def run_ska_allpairs(
    sample_records: List[Tuple[str, Path]],
    ska_root: Path,
    out_prefix: str,
    k: int = 31,
    threads: int | None = None,
):
    """
    Convenience wrapper:
      1) write file list
      2) ska build -> <out_prefix>.skf
      3) ska distance -> <out_prefix>.distances.tsv
    Returns (distances_tsv or None, combined_stdout).
    """
    filelist = ska_root / f"{out_prefix}.filelist.txt"
    write_filelist(sample_records, filelist)
    prefix = ska_root / out_prefix
    rc_b, out_b, skf = ska2_build(prefix, filelist, k=k)

    if rc_b != 0 or not skf.exists():
        return None, out_b

    out_tsv = ska_root / f"{out_prefix}.distances.tsv"
    rc_d, out_d = ska2_distance(skf, out_tsv, threads=threads)
    if rc_d != 0:
        return None, out_b + out_d
    return out_tsv, out_b + out_d
