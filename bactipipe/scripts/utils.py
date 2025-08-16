
import os
import re
import sys
import subprocess
import textwrap3
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from colorama import Fore, Style, init

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
    Computer Host Name: {hostname}
    Computer IP Address: {ip_address}
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
    
