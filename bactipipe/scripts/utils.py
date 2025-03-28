from datetime import datetime
from colorama import Fore, Style, init
import subprocess
import textwrap3
import sys
import os
from tqdm import tqdm

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
    
