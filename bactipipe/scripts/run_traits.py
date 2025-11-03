#!/usr/bin/env python3
import os, sys, argparse, json, subprocess, shlex
import multiprocessing as mp
from datetime import datetime
from pathlib import Path

# Reuse your existing BactiPipe utilities
from bactipipe.scripts.utils import time_print, logger as file_logger

from bactipipe.scripts import (
    traits_amr,
    traits_vf,
    traits_reconcile,
    traits_report,
    traits_db,
)

VIRAMR_ENV = os.getenv("VIRAMR_ENV", "viramr")
ABRICATE_ENV = os.getenv("ABRICATE_ENV", "abricate")

SUPPORTED_VF_ORGS = {
    "Escherichia",
    "Staphylococcus_aureus",
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Listeria",
}

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs['width'] = 100
        super().__init__(*args, **kwargs)
    def _get_help_string(self, action):
        if action.default is not None and action.default != argparse.SUPPRESS:
            return super()._get_help_string(action)
        return action.help

def _conda_run(env_name: str, cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-lc", f"conda run -n {shlex.quote(env_name)} {cmd}"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

def _list_amrfinder_organisms() -> list[str]:
    """
    Return AMRFinderPlus organism tokens (+ 'Listeria'), parsed only from the
    'Available --organism options:' line, ignoring Software/Database/Running/timing lines.
    """
    cp = _conda_run(VIRAMR_ENV, "amrfinder --list_organisms || true")
    raw = (cp.stdout or "") + "\n" + (cp.stderr or "")

    orgs: set[str] = set()
    target_prefix = "Available --organism options:"
    captured_line = None
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(target_prefix):
            captured_line = line
            break

    if captured_line:
        payload = captured_line.split(":", 1)[-1]
        for tok in payload.split(","):
            tok = tok.strip()
            if tok and all(ch.isalpha() or ch == "_" for ch in tok):
                orgs.add(tok)
    else:
        # Defensive fallback in case output format changes
        for line in raw.splitlines():
            line = line.strip()
            if not line or "Possible organisms" in line:
                line = line.split("Possible organisms:", 1)[-1].strip() if "Possible organisms" in line else ""
            if not line:
                continue
            for piece in line.replace(",", " ").split():
                piece = piece.strip()
                if not piece:
                    continue
                if piece.lower() in {"amrfinder", "took", "seconds", "to", "complete",
                                     "available", "database", "software", "running:", "running"}:
                    continue
                if all(ch.isalpha() or ch == "_" for ch in piece):
                    orgs.add(piece)

    orgs.add("Listeria")
    return sorted(orgs)

def _format_organism_list(orgs: list[str]) -> str:
    """
    Format organism list with '*' marking VirulenceFinder supported taxa.
    """
    lines = []
    lines.append("# '*' indicates VirulenceFinder-supported organisms.")
    for o in orgs:
        mark = "*" if o in SUPPORTED_VF_ORGS else " "
        lines.append(f"{mark} {o}")
    return "\n".join(lines)

def parse_args():
    p = argparse.ArgumentParser(
        prog="bactipipe traits-detect",
        description="Detect AMR determinants and virulence genes from assembled genomes.",
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    req = p.add_argument_group("Required arguments")
    opt = p.add_argument_group("Optional arguments")

    # Make these optional at parse-time; we'll validate unless --organism-list is used.
    req.add_argument("--analyst",
        help="Name/initials of the analyst (required).")
    req.add_argument("--genomes_dir",
        help="Directory of *.fasta assemblies (required).")
    req.add_argument("--sample_sheet",
        help="TSV (no header): isolate\\tsample\\tspecimen. Only first two columns used. "
             "(required).")
    req.add_argument("--run-name",
        help="Run name identifier, e.g., 250101_10_SEQ10_XYZ (required).")

    opt.add_argument("--accession", help="Run accession identifier (e.g., P25XXXXX)", default="Unknown")
    opt.add_argument("--outdir", help="Output directory. Default: ./Detect_<run-name>/<accession>")

    # Database root shared by CGE tools; expects subdirs like virulencefinder_db, resfinder_db, etc.
    opt.add_argument("--db-dir",
                     default=os.getenv("BACTIPIPE_DB_DIR"),
                     help="Root directory that contains virulencefinder_db, resfinder_db, etc.")
    req.add_argument("--organism",
                     help="Organism name for tool selection and AMRFinderPlus (-O) (required). Check spelling with --organism-list. if not listed, use genus name only.")
    opt.add_argument("--organism-list",
                     action="store_true",
                     help="Print organisms supported by AMRFinderPlus & VirulenceFinder and exit.")

    opt.add_argument("--threads", type=int, default=8, help="CPU threads")
    # Escape % in help strings for argparse
    opt.add_argument("--min_identity", type=float, default=90.0, help="Minimum %% identity")
    opt.add_argument("--min_coverage", type=float, default=60.0, help="Minimum %% coverage")

    mode = opt.add_mutually_exclusive_group()
    mode.add_argument("--amr-only", action="store_true", help="Run AMR only")
    mode.add_argument("--virulence-only", action="store_true", help="Run Virulence only")

    # Tool selectors (AMR) — AMRFinder runs by default; ABRicate toggles add-on if flagged
    opt.add_argument("--amrfinder", action="store_true",
                     help="(Optional) Explicitly enable AMRFinderPlus (on by default unless --virulence-only).")
    opt.add_argument("--abricate-resfinder", action="store_true",
                     help="Enable ABRicate/ResFinder (in addition to AMRFinder).")
    opt.add_argument("--abricate-card", action="store_true",
                     help="Enable ABRicate/CARD (in addition to AMRFinder).")

    # Tool selectors (Virulence) — if neither is given, choose based on organism
    opt.add_argument("--virulencefinder", action="store_true",
                     help="Enable VirulenceFinder (DTU Python tool).")
    opt.add_argument("--abricate-vfdb", action="store_true",
                     help="Enable ABRicate/VFDB (Set A). If provided, will be used regardless of the organism.")

    opt.add_argument("-h","--help", action="help", help="Show this help message and exit")
    return p.parse_args()

def _read_sample_sheet(sample_sheet: str, genomes_dir: str) -> list[tuple[str, str, str]]:
    """
    Returns a list of (isolate_id, sample_display, fasta_path).
    - sample_sheet columns: isolate \t sample \t specimen  (no header)
    - fasta_path = genomes_dir/{isolate}.fasta
    Only rows where fasta exists are kept.
    """
    rows: list[tuple[str, str, str]] = []
    missing: list[str] = []
    with open(sample_sheet, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            isolate = parts[0].strip()
            sample = parts[1].strip()
            fasta = os.path.join(genomes_dir, f"{isolate}.fasta")
            if os.path.isfile(fasta):
                rows.append((isolate, sample, fasta))
            else:
                missing.append(isolate)
    if missing:
        sys.stderr.write(f"WARNING: FASTA not found for isolates (skipped): {', '.join(missing)}\n")
    if not rows:
        sys.stderr.write("ERROR: No valid rows found in --sample_sheet (no matching FASTA files).\n")
        sys.exit(2)
    return rows

def _run_one(args_pack):
    isolate_id, sample_display, fasta, cfg = args_pack
    out_root = cfg["outdir"]
    log = os.path.join(out_root, "logs", f"{isolate_id}.log")
    os.makedirs(os.path.dirname(log), exist_ok=True)

    def logfn(msg, level="Norm", mode="timestamp"):
        file_logger(log, msg, message_type=level, mode=mode)

    logfn(f"Starting traits-detect for isolate={isolate_id} (report name={sample_display})")
    os.makedirs(os.path.join(out_root, "raw"), exist_ok=True)

    amr_rows = []
    vf_rows = []

    # Internally, we name outputs by isolate_id (matches FASTA and raw/ filenames)
    if not cfg["vir_only"]:
        amr_rows = traits_amr.run_amr_for_sample(
            sample=isolate_id,
            fasta=fasta,
            out_root=out_root,
            min_id=cfg["min_id"],
            min_cov=cfg["min_cov"],
            enable_amrfinder=cfg["amrfinder"],
            enable_resfinder=cfg["resfinder"],
            enable_card=cfg["card"],
            threads=cfg["threads"],
            organism=cfg["organism"],
            logger_fn=logfn,
        )

    if not cfg["amr_only"]:
        vf_rows = traits_vf.run_vf_for_sample(
            sample=isolate_id,
            fasta=fasta,
            out_root=out_root,
            min_id=cfg["min_id"],
            min_cov=cfg["min_cov"],
            enable_vfinder=cfg["virulencefinder"],  # DTU Python tool
            enable_vfdb=cfg["vfdb"],                # ABRicate VFDB
            threads=cfg["threads"],
            logger_fn=logfn,
            db_root=cfg["db_root"],
        )

    merged = traits_reconcile.merge_for_sample(isolate_id, amr_rows, vf_rows)

    # Per-sample TSVs are named by isolate_id to match raw outputs
    traits_report.write_sample_tsvs(isolate_id, merged, out_root)

    # Return using the display name as the key for reporting
    return sample_display, merged

def main():
    args = parse_args()

    # --- short-circuit: just list organisms and exit ---
    if args.organism_list:
        orgs = _list_amrfinder_organisms()
        print(_format_organism_list(orgs))
        sys.exit(0)

    # --- validate required args only for real runs ---
    missing = []
    if not args.analyst:      missing.append("--analyst")
    if not args.genomes_dir:  missing.append("--genomes_dir")
    if not args.sample_sheet: missing.append("--sample_sheet")
    if not args.organism:    missing.append("--organism")
    if not args.run_name:     missing.append("--run-name")

    if missing:
        sys.stderr.write("ERROR: missing required arguments: " + ", ".join(missing) + "\n")
        sys.exit(2)

    if not args.outdir:
        args.outdir = os.path.join(f"Detect_{args.run_name}", args.accession)

    os.makedirs(args.outdir, exist_ok=True)

    # --- discover samples from sample_sheet (maps isolate -> sample display) ---
    sample_rows = _read_sample_sheet(args.sample_sheet, args.genomes_dir)
    # sample_rows: list of (isolate_id, sample_display, fasta_path)

    # --- decide default tool selection per your logic ---
    # AMR defaults: AMRFinderPlus ON (unless virulence-only); ABRicate toggles add-on if flagged
    amr_amrfinder = (not args.virulence_only)  # default ON when AMR track runs
    if args.amrfinder:
        amr_amrfinder = True  # explicit still true

    amr_resfinder = bool(args.abricate_resfinder)
    amr_card = bool(args.abricate_card)

    # Virulence defaults:
    vf_vfinder = bool(args.virulencefinder)
    vf_vfdb = bool(args.abricate_vfdb)
    if not vf_vfinder and not vf_vfdb:
        if args.organism in SUPPORTED_VF_ORGS:
            vf_vfinder = True
        else:
            vf_vfdb = True

    # --- versions / db tracking ---
    versions = traits_db.ensure_and_collect_versions(
        want_amr=not args.virulence_only,
        want_vf=not args.amr_only,
        want_resfinder=amr_resfinder,
        want_card=amr_card,
        want_vfinder=vf_vfinder,
        want_vfdb=vf_vfdb,
        db_root=args.db_dir,
    )
    versions["thresholds"] = {
        "min_identity": float(args.min_identity),
        "min_coverage": float(args.min_coverage),
    }
    with open(os.path.join(args.outdir, "versions.json"), "w") as f:
        json.dump(versions, f, indent=2)

    cfg = dict(
        outdir=os.path.abspath(args.outdir),
        db_root=args.db_dir,
        min_id=float(args.min_identity),
        min_cov=float(args.min_coverage),
        threads=int(args.threads),

        amr_only=bool(args.amr_only),
        vir_only=bool(args.virulence_only),

        # Final tool toggles after default logic
        amrfinder=bool(amr_amrfinder),
        resfinder=bool(amr_resfinder),
        card=bool(amr_card),

        virulencefinder=bool(vf_vfinder),
        vfdb=bool(vf_vfdb),

        organism=args.organism,
        analyst=args.analyst,
        run_name=args.run_name,
        accession=args.accession,
    )

    start = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_print(f"traits-detect started at {start}", "Header")
    os.makedirs(os.path.join(cfg["outdir"], "raw"), exist_ok=True)
    os.makedirs(os.path.join(cfg["outdir"], "reports"), exist_ok=True)

    # Multiprocessing over *only* those isolates specified in the sample sheet
    pool = mp.Pool(processes=min(len(sample_rows), cfg["threads"]))
    try:
        results = pool.map(_run_one, [(iso, samp, fa, cfg) for (iso, samp, fa) in sample_rows])
    finally:
        pool.close(); pool.join()

    # --- unified run-level reporting ---
    # Key reports by the 'sample' (display name from column 2)
    all_merged = {sample_display: merged for (sample_display, merged) in results}

    # Unified TSV for convenience
    traits_report.write_run_summary(all_merged, cfg["outdir"])

    # Single PDF for the entire run
    traits_report.render_run_pdf(
        out_root=cfg["outdir"],
        run_name=cfg["run_name"],
        accession=cfg["accession"],
        tech_name=cfg["analyst"],
        versions=versions,
        thresholds=versions.get("thresholds", {}),
        all_merged=all_merged,
        pipeline_title="Virulence and AMR Gene Detection",
        run_amr=(not args.virulence_only),
        run_vf=(not args.amr_only),
    )
    # Also emit a consolidated matrix TSV (no header)
    accession = cfg["accession"]
    traits_report.write_consolidated_tsv(all_merged, cfg["outdir"], f"{accession}_vir_amr_detect.tsv")


    time_print("traits-detect completed.", "Header")

if __name__ == "__main__":
    main()
