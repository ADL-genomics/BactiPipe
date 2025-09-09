
import os
import re
import sys
import boto3
import shutil
import logging
import textwrap3
import subprocess
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from colorama import Fore, Style, init
from typing import Iterable, Callable, Optional, Dict, List, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

import io, gzip, threading, math
from Bio import SeqIO


_LEVEL_MAP = {
    "Norm": logging.INFO,
    "Pass": logging.INFO,
    "Header": logging.INFO,
    "Warn": logging.WARNING,
    "Fail": logging.ERROR,
}

def stringwraper(text, width, s_type):
    lines = textwrap3.wrap(text, width=width, break_long_words=True)
    new_lines = []
    if s_type == "plain":
        new_lines = lines
    elif s_type == "command":
        for line in lines[:-1]:
            new_lines.append(line + "\\")
        new_lines.append(lines[-1])
    return new_lines

def time_logger(logfile, message, message_type="", s_type='plain'):
    # logfile is an open file object
    # Adds a timestamp to the message and writes it to the logfile
    tstamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    # message = f"{tstamp} - {message}"
    messag_len = len(tstamp) + 2
    indent = " " * messag_len
    wrap_width = 70#70 - messag_len

    lines = stringwraper(message, width=wrap_width, s_type=s_type)
    if message_type == 'Header':
        logfile.write('\n')
        logfile.write(f'{tstamp} - {message}\n')
        logfile.write("="*len(f'{tstamp} - {message}') + '\n\n')
    elif message_type in ['Fail', 'Pass', 'Norm']:
        logfile.write(f'{tstamp} - {lines[0]}\n')
        for line in lines[1:]:
            logfile.write(indent + line + '\n')
    else:
        logfile.write(message + '\n')

def simple_logger(logfile, message, message_type="", s_type='plain'):
    # This allows for the logging without a timestamp
    if message_type == 'Header':
        logfile.write('\n')
        logfile.write(message + '\n')
        logfile.write("="*len(message) + '\n\n')
    else:
    	logfile.write(message + '\n')

def logger(logfile, message="", message_type="Norm", mode="timestamp", s_type='plain'):
    # Mode options: timestamp, simple
    logfile = os.path.abspath(logfile)
    if not os.path.exists(os.path.dirname(logfile)):
        os.makedirs(os.path.dirname(logfile))
    if mode == "timestamp":
        with open(logfile, 'a') as log:
            time_logger(logfile=log, message=message, message_type=message_type, s_type=s_type)
    elif mode == "simple":
        with open(logfile, 'a') as log:
            simple_logger(logfile=log, message=message, message_type=message_type, s_type=s_type)


def time_print(message, message_type="", s_type='plain'):
    tstamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    # message = f"{tstamp} - {message}"
    messag_len = len(tstamp) + 2
    indent = " " * messag_len
    wrap_width = 140 - messag_len

    lines = stringwraper(message, width=wrap_width, s_type=s_type)
  
    if message_type == 'Fail':
        print(Fore.RED + f'{tstamp} - {lines[0]}' + Style.RESET_ALL)  
        for line in lines[1:]:
            print(Fore.RED + indent + line + Style.RESET_ALL)
        # print('\n')

    elif message_type == 'Pass':
        print(Fore.GREEN + f'{tstamp} - {lines[0]}' + Style.RESET_ALL)  
        for line in lines[1:]:
            print(Fore.GREEN + indent + line + Style.RESET_ALL)
        # print('\n')

    elif message_type == 'Header':
        print('')
        print("\033[1m" + Fore.YELLOW + f'{tstamp} - {message}' + Style.RESET_ALL)  
        print("\033[1m" + Fore.YELLOW + "="*len(f'{tstamp} - {message}') + '\n' + Style.RESET_ALL)  
    else:
        print(f'{tstamp} - {lines[0]}')
        for line in lines[1:]:
            print(indent + line )
            # print('\n')

def simple_print(message, message_type="", s_type='plain'):
    if message_type == 'Fail':
        print(Fore.RED + message.lstrip() + Style.RESET_ALL)  

    elif message_type == 'Pass':
        print(Fore.GREEN + message.lstrip() + Style.RESET_ALL)  

    elif message_type == 'Header':
        print('\n')
        print("\033[1m" + Fore.YELLOW + message.lstrip() + Style.RESET_ALL)  # Green for positive
        print("\033[1m" + Fore.YELLOW + "="*len(message.lstrip()) + '\n' + Style.RESET_ALL)  # Green for positive
    else:
        print(message)  

        
def check_required(programs, python_packages):
    missing_programs = []
    missing_packages = []
    
    for program in programs:
        result = subprocess.run(['which', program], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            missing_programs.append(program)
    
    for package in python_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_programs or missing_packages:
        if missing_programs:
            print("The following required programs are not installed:")
            for program in missing_programs:
                print(f"- {program}")
        if missing_packages:
            print("The following required Python packages are not installed:")
            for package in missing_packages:
                print(f"- {package}")
        print("Please install the missing programs and packages and try again.")
        sys.exit(1)

def pipeheader(date, tech_name, hostname, ip_address, run_name, sample_list, raw_reads, outDir, cpus):
    header = f'''
    =======================================
    Penn State Animal Diagnostic Laboratory
    Bacteria WGS QC Pipeline v.1.0.0
    Run name: {run_name}
    Date: {date}
    Analysis by: {tech_name}
    Server Host Name: {hostname}
    =======================================
    '''
    header = header.split('\n')
     
    run_info = f'''
    Run Parameters:
    ===============
    \tRun name: {run_name}
    \tSample sheet: {sample_list}
    \tInput directory (fastq reads): {raw_reads}
    \tOutput directory (results): {os.path.abspath(outDir)}
    \tNumber of threads: {cpus}
    '''
    run_info = run_info.split('\n')
    return header, run_info

def excel_reader(filepath: str, platform: str = "nanopore"):
    # — normalize function treats NaN as empty and collapses whitespace/linebreaks —
    def normalize(x):
        if pd.isna(x):
            return ""
        s = str(x).strip().strip('"')
        return re.sub(r'\s+', ' ', s)

    plat = platform.lower()
    if plat == "nanopore":
        needed = [
            "Accession Number",
            "Host and Source",
            "Bacteria Species",
            "gDNA Prep ID",
            "Native Barcode Index",
        ]
        idx_col = "Native Barcode Index"
        pattern = re.compile(r"^(?:NB\d{2}|barcode\d{2})$", re.IGNORECASE)
        header_lines = 1

    elif plat == "illumina":
        needed = [
            "Well",
            "Accession Number",
            "Host and Source",
            "Bacteria Species",
            "gDNA Prep ID",
            "UD Set A Index Well",
        ]
        idx_col = "Well"
        pattern = re.compile(r"^[A-H](?:0[1-9]|1[0-2])$")
        header_lines = 2

    else:
        raise ValueError("`platform` must be 'nanopore' or 'illumina'")

    # 1) Read every cell as raw (no header)
    raw = pd.read_excel(filepath, sheet_name=0, header=None, engine="openpyxl")
    n_rows, n_cols = raw.shape

    # 2) Find the header block
    header_row = None
    header_vals = None

    for i in range(n_rows - header_lines + 1):
        if header_lines == 1:
            row_vals = [normalize(x) for x in raw.iloc[i]]
        else:
            top = [normalize(x) for x in raw.iloc[i]]
            bot = [normalize(x) for x in raw.iloc[i + 1]]
            row_vals = [top[j] or bot[j] for j in range(n_cols)]
        if all(col in row_vals for col in needed):
            header_row = i
            header_vals = row_vals
            break

    if header_row is None:
        raise ValueError(f"Could not locate the header block for platform '{platform}'.")

    # 3) Now safely compute where data starts
    data_start = header_row + header_lines

    # 4) Map each needed name to its column index
    col_positions = [header_vals.index(col) for col in needed]

    # 5) Fix merged‐cell misalignment:
    #    if the “header” column j is blank for all data rows
    #    but j+1 has real values, shift to j+1
    for idx, j in enumerate(col_positions):
        block = raw.iloc[data_start : data_start + 10, j]
        # is it entirely blank?
        if all(pd.isna(x) or str(x).strip() == "" for x in block):
            # does the next column have any data?
            if j + 1 < n_cols and any(
                not (pd.isna(x) or str(x).strip() == "")
                for x in raw.iloc[data_start : data_start + 10, j + 1]
            ):
                col_positions[idx] = j + 1

    # 6) Slice out exactly those columns, build a DataFrame
    data_block = raw.iloc[data_start:, col_positions]
    df = pd.DataFrame(data_block.values, columns=needed)

    # 7) Clean up every cell
    df = df.apply(lambda col: col.map(normalize))

    # 8) Keep only rows where the index‐column matches the regex
    df = df[df[idx_col].astype(bool) & df[idx_col].str.match(pattern)]

    # 9) Drop any row missing the two critical fields
    df = df[df["Accession Number"].astype(bool) & df["Bacteria Species"].astype(bool)]

    # 10) Return a list of TAB‐joined strings

    if plat == "nanopore":
        out = df.apply(
            lambda row: f"{row['Accession Number']}\t{row['Bacteria Species']}\t{row[idx_col]}",
            axis=1
        ).tolist()
    elif plat == "illumina":
        out = df.apply(
            lambda row: f"{row['Accession Number']}\t{row['Bacteria Species']}",
            axis=1
        ).tolist()
    
    return out

def download_s3(bucket: str, prefix: str, target_root: str) -> str:
    s3 = boto3.client("s3")
    norm_prefix = prefix.strip("/")
    base_name = norm_prefix.split("/")[-1]
    local_base = os.path.join(target_root, base_name)
    os.makedirs(local_base, exist_ok=True)

    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=norm_prefix + "/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            rel_path = key[len(norm_prefix):].lstrip("/")
            if not rel_path:
                continue
            dest_path = os.path.join(local_base, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            s3.download_file(bucket, key, dest_path)

    return local_base


init(autoreset=True)

class TqdmMinutes(tqdm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._rate_format = "{rate:.2f} min/genome"

    def format_dict(self):
        d = super().format_dict
        rate_min = (self.format_dict['elapsed'] / self.n) / 60 if self.n else 0
        d['rate'] = rate_min
        return d

    def format_meter(self, n, total, elapsed, ncols=None):
        rate_min = (elapsed / n) / 60 if n else 0
        return super().format_meter(n, total, elapsed, ncols).replace("it/s", f"{rate_min:.2f} min/genome")
    
# Compress files at the end of run
def compress_qc_fastqs(
    qc_out: str,
    genomes_dir: str,
    cpus_per_sample: int = 4,
    require_assembly: bool = True,
    recursive: bool = False,
    remove_source: bool = True,
    logger_fn: Optional[Callable[[str, str], None]] = None,
) -> Dict[str, int]:
    """
    Compress uncompressed FASTQ files under qc_out/<sample> to .gz in parallel (threads).

    - Uses pigz if available (multi-threaded), else gzip.
    - Atomic write: *.gz.tmp -> os.replace(...).
    - Verifies integrity (pigz/gzip -t) before deleting source.
    - Avoids CPU oversubscription: workers ~= total_cpus // threads_per_file.
    - Logs via standard `logging` so records go into your main log file.
      Optionally accepts `logger_fn(message, level)`; if provided, it’s used first.

    Returns: {"total": N, "ok": K, "failed": F, "skipped": S}
    """
    log = logging.getLogger(__name__)

    def _log(msg: str, level: str = "Norm") -> None:
        if logger_fn:
            try:
                logger_fn(msg, level)   # e.g., compat_logger(..., message_type=level)
                return
            except TypeError:
                try:
                    logger_fn(msg)      # accept simple print-like callables too
                    return
                except Exception:
                    pass
        log.log(_LEVEL_MAP.get(level, logging.INFO), msg)

    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        _log("Neither pigz nor gzip found on PATH; skipping compression.", "Warn")
        return {"total": 0, "ok": 0, "failed": 0, "skipped": 0}

    # ---- Gather jobs ----------------------
    jobs: List[Tuple[str, str]] = []  # (in_fastq, out_gz)
    skipped = 0

    for sample in sorted(os.listdir(qc_out)):
        sample_dir = os.path.join(qc_out, sample)
        if not os.path.isdir(sample_dir):
            continue

        if require_assembly:
            asm_path = os.path.join(genomes_dir, f"{sample}.fasta")
            if not (os.path.exists(asm_path) and os.path.getsize(asm_path) > 0):
                skipped += 1
                continue

        if recursive:
            for root, _, files in os.walk(sample_dir):
                for name in files:
                    if name.endswith(".fastq"):         # add ".fq" here if you want
                        in_path = os.path.join(root, name)
                        out_gz = in_path + ".gz"
                        if not os.path.exists(out_gz):
                            jobs.append((in_path, out_gz))
        else:
            for name in os.listdir(sample_dir):
                if name.endswith(".fastq"):
                    in_path = os.path.join(sample_dir, name)
                    out_gz = in_path + ".gz"
                    if not os.path.exists(out_gz):
                        jobs.append((in_path, out_gz))

    if not jobs:
        _log("No uncompressed FASTQ files found to compress.", "Norm")
        return {"total": 0, "ok": 0, "failed": 0, "skipped": skipped}

    # ---- Parallelism budgeting --------------------------------------------
    try:
        threads_per_file = max(1, int(cpus_per_sample))
    except Exception:
        threads_per_file = 4
    if os.path.basename(compressor) != "pigz":
        threads_per_file = 1  # gzip is single-threaded

    total_cpus = os.cpu_count() or 1
    max_workers = max(1, min(len(jobs), max(1, total_cpus // threads_per_file)))

    _log(
        f"Compressing {len(jobs)} FASTQ files with {os.path.basename(compressor)} "
        f"(workers={max_workers}, threads/file={threads_per_file}).",
        "Norm",
    )

    # ---- Helpers -----------------------------------------------------------
    def _integrity_ok(gz_path: str) -> bool:
        tester = shutil.which("pigz") or shutil.which("gzip")
        if not tester:
            return True
        r = subprocess.run(
            [tester, "-t", gz_path],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return r.returncode == 0

    def _compress_one(in_path: str, out_gz: str) -> str:
        gz_tmp = out_gz + ".tmp"
        # Remove stale tmp if present
        try:
            if os.path.exists(gz_tmp):
                os.remove(gz_tmp)
        except Exception:
            pass

        cmd = [compressor]
        if os.path.basename(compressor) == "pigz" and threads_per_file > 1:
            cmd += ["-p", str(threads_per_file)]
        cmd += ["-c", in_path]

        # Primary attempt
        try:
            with open(gz_tmp, "wb") as gz_out:
                subprocess.run(
                    cmd,
                    check=True,
                    stdin=subprocess.DEVNULL,
                    stdout=gz_out,
                    stderr=subprocess.DEVNULL,
                )
        except subprocess.CalledProcessError:
            # Fallback to gzip if pigz failed
            if os.path.basename(compressor) == "pigz":
                gz = shutil.which("gzip")
                if not gz:
                    raise
                with open(gz_tmp, "wb") as gz_out:
                    subprocess.run(
                        [gz, "-c", in_path],
                        check=True,
                        stdin=subprocess.DEVNULL,
                        stdout=gz_out,
                        stderr=subprocess.DEVNULL,
                    )
            else:
                raise

        # Verify gzip integrity
        if not _integrity_ok(gz_tmp):
            try:
                os.remove(gz_tmp)
            except FileNotFoundError:
                pass
            raise RuntimeError(f"Compression test failed for {in_path}")

        os.replace(gz_tmp, out_gz)  # atomic publish

        if remove_source:
            try:
                os.remove(in_path)
            except FileNotFoundError:
                pass

        return out_gz

    # ---- Run in thread pool -----------------------------------------------
    ok = 0
    failed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(_compress_one, src, dst): (src, dst) for (src, dst) in jobs}
        for i, fut in enumerate(as_completed(futs), start=1):
            src, _dst = futs[fut]
            try:
                fut.result()
                ok += 1
                if i % 10 == 0 or i == len(jobs):
                    _log(f"Compressed {ok}/{len(jobs)} files.", "Norm")
            except Exception as e:
                failed += 1
                _log(f"Compression failed for {src}: {e}", "Warn")

    if failed:
        _log(f"End-of-pipeline compression finished with {failed} failure(s).", "Warn")
    else:
        _log("End-of-pipeline compression finished successfully.", "Pass")

    return {"total": len(jobs), "ok": ok, "failed": failed, "skipped": skipped}

# ---------- fast, safe open for SeqIO (supports .gz via pigz/gzip -dc) ----------
def _open_text_fastq_for_seqio(path: str, pigz_threads: int = 1) -> Tuple[io.TextIOBase, Optional[subprocess.Popen]]:
    """
    Returns (text_handle, proc). If proc is not None, caller should wait() it after use.
    """
    proc = None
    if path.endswith(".gz"):
        dec = shutil.which("pigz") or shutil.which("gzip")
        if dec:
            # robust int for pigz threads
            try:
                t = max(1, int(pigz_threads))
            except Exception:
                t = 1
            cmd = [dec]
            if os.path.basename(dec) == "pigz" and t > 1:
                cmd += ["-p", str(t)]
            cmd += ["-dc", path]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            assert proc.stdout is not None
            fh = io.TextIOWrapper(io.BufferedReader(proc.stdout, 1024 * 1024),
                                  encoding="ascii", newline="")
            return fh, proc
        else:
            fh = io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", newline="")
            return fh, None
    else:
        return open(path, "rt", encoding="ascii", newline=""), None


# =========================== API 1: fastq_metrics ============================

def fastq_metrics(path: str, use_external_decompress: bool = True, pigz_threads=1):
    """
    Streaming metrics from FASTQ/FASTQ.GZ using Biopython/SeqIO.
    Returns a dict with:
      reads, total_bases,
      mean_q_read_np   -> NanoPlot/MinKNOW-style mean read quality (per-read prob avg -> Phred)
      mean_q_base_np   -> Phred of base-weighted mean error prob across all bases
      mean_q_per_base_arith -> arithmetic mean of Phred across bases (legacy/comparison)
      mean_q_per_read_arith -> arithmetic mean of per-read mean Phred (legacy/comparison)
    """
    if path.endswith(".gz") and use_external_decompress:
        fh, proc = _open_text_fastq_for_seqio(path, pigz_threads=pigz_threads)
    else:
        fh, proc = open(path, "rt", encoding="ascii", newline=""), None

    reads = 0
    total_bases = 0

    # NanoPlot/MinKNOW accumulators
    p_read_sum = 0.0           # sum over reads of mean per-base error prob
    p_base_sum = 0.0           # sum over all bases of error prob (i.e., Σ_read p_read * L)

    # Arithmetic means (for comparison with FastQC-like summaries)
    qsum_all = 0               # Σ over all bases of Phred
    read_mean_qsum = 0.0       # Σ over reads of (mean Phred per read)

    try:
        for rec in SeqIO.parse(fh, "fastq"):   # PHRED+33 integers
            quals = rec.letter_annotations.get("phred_quality")
            if not quals:
                continue
            L = len(quals)
            if L <= 0:
                continue

            # NanoPlot/MinKNOW-style per-read mean error prob
            p_read = sum(10 ** (-q / 10.0) for q in quals) / L
            p_read_sum += p_read
            p_base_sum += p_read * L

            # Arithmetic accumulators
            qsum = sum(quals)
            qsum_all += qsum
            read_mean_qsum += (qsum / L)

            reads += 1
            total_bases += L
    finally:
        fh.close()
        if proc is not None:
            try: proc.wait()
            except Exception: pass

    # Convert prob means back to Phred (guard against zeros)
    mean_q_read_np = (-10.0 * math.log10(p_read_sum / reads)) if reads and p_read_sum > 0 else 0.0
    mean_q_base_np = (-10.0 * math.log10(p_base_sum / total_bases)) if total_bases and p_base_sum > 0 else 0.0

    # Arithmetic means
    mean_q_per_base_arith = (qsum_all / total_bases) if total_bases else 0.0
    mean_q_per_read_arith = (read_mean_qsum / reads) if reads else 0.0

    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_q_read_np": mean_q_read_np,
        "mean_q_base_np": mean_q_base_np,
        "mean_q_per_base_arith": mean_q_per_base_arith,
        "mean_q_per_read_arith": mean_q_per_read_arith,
    }


# ================= API 2: run_filtlong_and_write_with_metrics ================

def filtlong_with_metrics(filtlong_cmd, out_fastq_path: str):
    """
    Run filtlong, stream its FASTQ output, write UNCOMPRESSED out_fastq_path,
    and compute NanoPlot/MinKNOW-style metrics in one pass.
    Returns the same keys as fastq_metrics(...).
    """
    p = subprocess.Popen(
        filtlong_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,       # drain stderr so it can't block
        stdin=subprocess.DEVNULL,
        bufsize=1024 * 1024,
        text=False,
    )
    assert p.stdout is not None

    # drain stderr in background
    def _drain_err(stream):
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break

    t = threading.Thread(target=_drain_err, args=(p.stderr,), daemon=True)
    t.start()

    # Wrap stdout for text parsing
    fh = io.TextIOWrapper(io.BufferedReader(p.stdout, 1024 * 1024),
                          encoding="ascii", newline="")

    reads = 0
    total_bases = 0
    p_read_sum = 0.0
    p_base_sum = 0.0
    qsum_all = 0
    read_mean_qsum = 0.0

    with open(out_fastq_path, "wt", encoding="ascii", newline="") as out:
        for rec in SeqIO.parse(fh, "fastq"):
            quals = rec.letter_annotations.get("phred_quality")
            if not quals:
                continue
            L = len(quals)
            if L <= 0:
                continue

            # NanoPlot-style accumulators
            p_read = sum(10 ** (-q / 10.0) for q in quals) / L
            p_read_sum += p_read
            p_base_sum += p_read * L

            # Arithmetic accumulators (optional, for comparison)
            qsum = sum(quals)
            qsum_all += qsum
            read_mean_qsum += (qsum / L)

            reads += 1
            total_bases += L

            # Write normalized 4-line FASTQ
            SeqIO.write(rec, out, "fastq")

    ret = p.wait()
    t.join()
    if ret != 0:
        raise subprocess.CalledProcessError(ret, filtlong_cmd)

    mean_q_read_np = (-10.0 * math.log10(p_read_sum / reads)) if reads and p_read_sum > 0 else 0.0
    mean_q_base_np = (-10.0 * math.log10(p_base_sum / total_bases)) if total_bases and p_base_sum > 0 else 0.0
    mean_q_per_base_arith = (qsum_all / total_bases) if total_bases else 0.0
    mean_q_per_read_arith = (read_mean_qsum / reads) if reads else 0.0

    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_q_read_np": mean_q_read_np,
        "mean_q_base_np": mean_q_base_np,
        "mean_q_per_base_arith": mean_q_per_base_arith,
        "mean_q_per_read_arith": mean_q_per_read_arith,
    }


## For merging files
def merge_gz_fastqs(raw_folder: str, output_gz: str, *, verify: bool = False) -> str:
    """
    Merge all *.fastq.gz in raw_folder into output_gz using a single `cat`.
    Writes to <output_gz>.tmp then atomically renames. Returns output path.
    """
    if not output_gz.endswith(".gz"):
        raise ValueError("output_gz must end with .gz")

    out_dir = os.path.dirname(os.path.abspath(output_gz)) or "."
    os.makedirs(out_dir, exist_ok=True)

    tmp_out = output_gz + ".tmp"

    # One fast cat. Important: command string is a single, quoted argument after -c.
    cmd = f'cat "{raw_folder}"/*.fastq.gz > "{tmp_out}"'
    subprocess.run(["/bin/sh", "-c", cmd], check=True)

    if verify:
        gz = shutil.which("pigz") or shutil.which("gzip")
        if gz:
            subprocess.run([gz, "-t", tmp_out], check=True)

    os.replace(tmp_out, output_gz)
    return output_gz

#S3 sources: merge_gz_fastqs_s3

def _parse_s3_uri(uri: str) -> Tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"Not an s3 URI: {uri}")
    without = uri[5:]
    bucket, _, key = without.partition("/")
    if not bucket:
        raise ValueError(f"Invalid s3 URI: {uri}")
    return bucket, key

def _list_s3_fastq_gz(session, bucket: str, prefix: str) -> List[Tuple[str, str]]:
    s3 = session.client("s3")
    paginator = s3.get_paginator("list_objects_v2")
    keys: List[Tuple[str, str]] = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            k = obj["Key"]
            if k.endswith((".fastq.gz", ".fq.gz")):
                keys.append((bucket, k))
    keys.sort(key=lambda bk: _natural_key(bk[1]))
    return keys

def merge_gz_fastqs_s3(
    s3_sources_or_prefix: Iterable[str] | str,
    output_path: str,
    aws_profile: Optional[str] = None,
    region: Optional[str] = None,
    verify: bool = False,
    logger_fn: Optional[Callable[[str, str], None]] = None,
) -> dict:
    """
    Fast merge of S3 .fastq.gz/.fq.gz objects -> local .gz or uncompressed FASTQ.
    Requires boto3.
    """
    def log(msg: str, level: str = "Norm"):
        if logger_fn:
            try:
                logger_fn(msg, level)
            except Exception:
                pass

    try:
        import boto3
    except ImportError:
        raise RuntimeError("boto3 is required for merge_gz_fastqs_s3.")

    session = boto3.Session(profile_name=aws_profile, region_name=region)

    # Build (bucket, key) list
    inputs: List[Tuple[str, str]] = []
    if isinstance(s3_sources_or_prefix, str):
        bucket, key = _parse_s3_uri(s3_sources_or_prefix)
        if key.endswith("/"):
            inputs = _list_s3_fastq_gz(session, bucket, key)
        else:
            inputs = [(bucket, key)]
    else:
        for uri in s3_sources_or_prefix:
            b, k = _parse_s3_uri(str(uri))
            inputs.append((b, k))
        inputs.sort(key=lambda bk: _natural_key(bk[1]))

    if not inputs:
        raise FileNotFoundError("No S3 .fastq.gz/.fq.gz objects found to merge.")

    s3 = session.client("s3")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    tmp_out = output_path + ".tmp"

    try:
        if output_path.endswith(".gz"):
            log(f"Merging {len(inputs)} S3 gz FASTQs -> {output_path} (raw concat).")
            total = 0
            with open(tmp_out, "wb") as out:
                for bucket, key in inputs:
                    obj = s3.get_object(Bucket=bucket, Key=key)
                    body = obj["Body"]
                    for chunk in iter(lambda: body.read(8 * 1024 * 1024), b""):
                        out.write(chunk)
                        total += len(chunk)
            if verify and not _gz_integrity_ok(tmp_out):
                raise RuntimeError("Gzip integrity test failed for merged output.")
            os.replace(tmp_out, output_path)
            return {"inputs": len(inputs), "bytes_written": total, "output": output_path}
        else:
            log(f"Merging {len(inputs)} S3 gz FASTQs -> {output_path} (decompress).")
            total = 0
            import gzip  # python gzip to decompress stream
            with open(tmp_out, "wb") as out:
                for bucket, key in inputs:
                    obj = s3.get_object(Bucket=bucket, Key=key)
                    body = obj["Body"]
                    # Wrap the streaming body in a file-like for gzip
                    with gzip.GzipFile(fileobj=body, mode="rb") as gz:
                        shutil.copyfileobj(gz, out, length=8 * 1024 * 1024)
            total = os.path.getsize(tmp_out)
            os.replace(tmp_out, output_path)
            return {"inputs": len(inputs), "bytes_written": total, "output": output_path}
    finally:
        try:
            if os.path.exists(tmp_out):
                os.remove(tmp_out)
        except Exception:
            pass
