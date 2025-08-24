import os
import sys
import argparse
import csv
import time
import socket
import shutil
import subprocess
from datetime import datetime
from bactipipe.scripts import find_organism
from bactipipe.scripts import process_data
from bactipipe.scripts.qualityProc import ProcQuality
from bactipipe.scripts.utils import time_print, simple_print, logger, pipeheader, excel_reader, download_s3
from importlib.resources import files as resource_files


# Argument Parsing
class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs['width'] = 100  # Set the width of the help message
        super().__init__(*args, **kwargs)

    def _get_help_string(self, action):
        if action.default is not None and action.default != argparse.SUPPRESS:
            return super()._get_help_string(action)
        return action.help

parser = argparse.ArgumentParser(
    description="Perform quality control on nanopore reads, assemble genomes for reads that pass QC, and check quality for assembled genomes.",
    formatter_class=CustomHelpFormatter,
    add_help=False  # Disable default help to add it manually
)

required_args = parser.add_argument_group('Required arguments')
optional_args = parser.add_argument_group('Optional arguments')


required_args.add_argument("-n", "--name",
                    help='Name of the person running the pipeline. Needs quotations to handle spaces.\
                        Example: -n "Jane Doe"', required=True)

required_args.add_argument("-r", "--run_name",
                    help="Run name (Run ID)", required=True)

required_args.add_argument("-d", "--fastq_dir",
                    help="Folder containing subfloders of raw fastq files from Illumina instrument.", required=True)

required_args.add_argument("-s", "--storage_type", choices=["local", "s3", "network_mount"], default="local",
                    help="Type of storage where the raw data is located. Default: local")

required_args.add_argument("-l", "--sample_sheet",
                    help="A two-column file with sample (column 1) and organism species (column2).", required=True)

optional_args.add_argument("-o", "--outdir",
                    help="Path to output directory. Default: current directory")

optional_args.add_argument ("-a", "--assembler", help="Genome assembler to use. Default: Spades",  choices=["Spades", "Unicycler"], default="Spades")

optional_args.add_argument("-t", "--threads",
                    help="Number of threads. Default: all available")



optional_args.add_argument("-h", "--help",
                           action="help",
                           help="Show this help message and exit.")

args = parser.parse_args()

sample_list = os.path.abspath(args.sample_sheet)
run_name = args.run_name
raw_reads = os.path.abspath(args.fastq_dir)
raw_reads = os.path.join(raw_reads, run_name)
tech_name = args.name

outDir = os.path.join(args.outdir, run_name)

def get_fastqs(sample, raw_reads):
    fastqs = []
    for item in os.listdir(raw_reads):
        if item.startswith(sample):
            for file in os.listdir(os.path.join(raw_reads, item)):
                if "_R1" in file and (file.endswith("fastq") or file.endswith("fastq.gz")):
                    fastq1 = os.path.join(raw_reads, item, file)
                elif "_R2" in file and (file.endswith("fastq") or file.endswith("fastq.gz")):
                    fastq2 = os.path.join(raw_reads, item, file)
    if not fastq2:  # If the data is single end
        fastqs.append(fastq1)
    else:
        fastqs.append(fastq1)
        fastqs.append(fastq2)
    return fastqs

# Run the pipeline

if not os.path.exists(outDir):
    os.makedirs(outDir)

if args.threads:
    cpus = args.threads
else:
    cpus = os.cpu_count()

static = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

log_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log = os.path.join(outDir, f"{log_time}~{run_name}_qc.log")

print("")

date = datetime.now().strftime("%Y-%m-%d")

hostname = socket.gethostname()
try:
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't even have to be reachable
        s.connect(('10.254.254.254', 1))
        ip_address = s.getsockname()[0]
    except Exception as e:
        ip_address = f"Unable to resolve hostname: {e}"
    finally:
        s.close()
except Exception as e:
    ip_address = f"Unable to resolve hostname: {e}"
    ip_address = "Unable to resolve hostname"

header_info = pipeheader(date, tech_name, hostname, ip_address, run_name, sample_list, raw_reads, outDir, cpus)

header = header_info[0]
run_info = header_info[1]
for line in header:
    simple_print(line)
    logger(log, line, mode="simple")
for line in run_info:
    simple_print(line)
    logger(log, line, mode="simple")

# Ensure that the sample list if formatted correctly
if not os.path.exists(sample_list):
    print(f"Sample list file {sample_list} does not exist.")
    logger(log, f"Sample list file {sample_list} does not exist.")
    sys.exit(1)
else:
    if sample_list.endswith('.xlsx') or sample_list.endswith('.xls'):
        try:
            sample_info = excel_reader(sample_list, "illumina")
        except Exception as e:
            print(f"Error reading Excel file {sample_list}: {e}")
            logger(log, f"Error reading Excel file {sample_list}: {e}")
            sys.exit(1)
    else:
        with open(sample_list, 'r') as sampL:
            sample_info = sampL.readlines()
            first_line = sample_info[0].strip().split('\t')
            if len(first_line) != 2:
                print("Sample list file must have two columns: sample, organism.")
                logger(log, "Sample list file must have two columns: sample, organism.")
                sys.exit(1)

# Load pathogenic bacteria information
bacteria = {}
org_list = resource_files('bactipipe.data').joinpath('pathogenic_bacteria.txt')

with open(org_list, 'r') as orgL:
    for line in orgL:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            name, size, taxo = parts[-1], parts[1], parts[0]
            bacteria[name] = [int(size), taxo]

# Validate the sample sheet
bad_organisms = []
bad_samples = []
# with open(sample_list, 'r') as sampL:
for line in sample_info:
    if line.startswith("#"):
        continue  # Skip comment lines
    sample, organism = line.strip().split('\t')
    logger(log, f"Validating sample: {sample}, organism: {organism}")
    if organism.strip() not in bacteria and organism.lower() not in ["organism", "unknown"]:
        time_print(f"WARNING: Organism [{organism.strip()}] for sample {sample} is not a valid organism name")
        logger(log, f"WARNING: Organism [{organism.strip()}] for sample {sample} is not a valid organism name")
        bad_organisms.append(f"{sample}:{organism}")
    else:
        print(f"Organism [{organism.strip()}] for sample {sample} is valid.")
    if not sample in ','.join(os.listdir(raw_reads)) and sample != "sample_ID":
        time_print(f"WARNING: Sample [{sample}] does not have corresponding fastq files", "Fail")
        logger(log, f"WARNING: Sample [{sample}] does not have corresponding fastq files")
        bad_samples.append(sample)

if bad_organisms:
    time_print("WARNING: Invalid organism names may cause QC failures.", "Fail")
    logger(log, "WARNING: Invalid organism names may cause QC failures.")

if bad_samples:
    b_s = ", ".join(bad_samples)
    time_print(f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.", "Fail")
    logger(log, f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.")

start_time = time.time()
logger(log, "PIPELINE STARTED", "Header")
time_print("PIPELINE STARTED", "Header")

# Count total number of samples
sample_number = sum(1 for line in sample_info if not line.startswith("#")) - len(bad_samples)

time_print(f"Total number of samples to be processed: {sample_number}\n")
logger(log, f"Total number of samples to be processed: {sample_number}\n")

qc_out = os.path.join(outDir, "qc_out")
if not os.path.exists(qc_out):
    os.makedirs(qc_out)

temp_dir = os.path.join(outDir, "temp")
if args.storage_type != "local":
    os.makedirs(temp_dir, exist_ok=True)
    logger(log, f"Storage @ {args.storage_type} - Created temporary directory for raw reads: {temp_dir}")
   
if args.storage_type == "s3":
    if not os.path.exists(temp_dir) or not os.listdir(temp_dir):
        time_print("Downloading raw reads from S3 bucket...")
        logger(log, "Downloading raw reads from S3 bucket...")
        # Download the files
        raw_reads = download_s3(bucket="adl-pathogen-bucket", prefix=f"ngslab/DATA/illuminaData/{run_name}", target_root=temp_dir)
    else:
        time_print("WARNING: Using existing temporary directory for raw reads.")
        logger(log, "WARNING:Using existing temporary directory for raw reads.")
        raw_reads = os.path.join(temp_dir, run_name)

elif args.storage_type == "network_mount":
    cp_cmd = f"cp -r {raw_reads} {temp_dir}"
    subprocess.run(cp_cmd, shell=True, check=True)
    raw_reads = os.path.join(temp_dir, os.path.basename(raw_reads))

#Write the temporary summary file that will be updated after CheckM analysis
temp_qc_summary = os.path.join(qc_out, "temp_qc_summary.tsv")
with(open(temp_qc_summary , 'w')) as qc_sum:
    writer = csv.writer(qc_sum, dialect='excel-tab')
    writer.writerow(["Sample",  "Mean_quality", "qc_verdict", "Expected organism", "Identified organism", "% Match", "Coverage Depth", "min_cov", "cov_verdict", "tax_confirm"])

    # with open(sample_list, 'r') as sampL:
    for line in sample_info:
        sample, organism = line.strip().split('\t')
        fastqs = get_fastqs(sample, raw_reads)
        if organism not in bacteria:
            sys_organism = "unknown"
            genome_size = None
            taxID = "N/A"
        else:
            sys_organism = organism
            genome_size = bacteria.get(organism, None)[0]
            taxID = bacteria.get(organism, None)[1]
        qc_file = os.path.join(qc_out, f"{sample}_metrics.txt")
        if not os.path.exists(qc_file) or os.path.getsize(qc_file) == 0:
            qualdata = ProcQuality(fastqs, genome_size, os.path.join(qc_out, sample), sample, logfile=log)

            qualdata.plot_quality_distribution()
            qualdata.plot_quality_metrics()

            avqc = qualdata.average_quality
            coverage = qualdata.coverage
        else:
            with open(qc_file, 'r') as qcFile:
                lines = qcFile.readlines()
            qual_line = next((l for l in lines if "Avg Quality:" in l), None)
            avqc = float(qual_line.split(":")[1].strip()) if qual_line else 0.0

            bases_line = next((l for l in lines if "Total bases:" in l), None)
            if bases_line is None:
                bases_line = next((l for l in lines if "Tot bases:" in l), None)

            total_bases = int(float(bases_line.split(":")[1].strip().replace(",", ""))) if bases_line else 0
            coverage = total_bases / genome_size if genome_size and total_bases > 0 else "N/A"

        taxonomy = [taxID, organism]

        if avqc >= 28:
            qc_verdict = 'Pass'
        else:
            qc_verdict = "Fail"

        # Coverage
        if avqc >= 30:
            min_cov = 30
        elif avqc >= 29:
            min_cov = 40
        elif avqc >= 28:
            min_cov = 50

        if organism == "Escherichia coli":
            min_cov += 10

        if coverage == "N/A":
            cov_verdict = "Pass"
        elif int(coverage) >= min_cov:
            cov_verdict = "Pass"
        else:
            cov_verdict = "Fail"

        simple_print(f'\t---> Mean quality score: {qc_verdict} --- {avqc:.2f}', qc_verdict)
        logger(log, f'\t---> Mean quality score: {qc_verdict} --- {avqc:.2f}', qc_verdict, mode="simple")

        if coverage != 'N/A':
            simple_print(f'\t---> Genome coverage: {cov_verdict} --- {coverage:.2f}X', cov_verdict)
            logger(log, f'\t---> Genome coverage: {cov_verdict} --- {coverage:.2f}X', cov_verdict, mode="simple")
            cov_display = f"{coverage:.2f}X"
        else:
            cov_display = "N/A"
            time_print(f'\t---> Genome coverage: {coverage}')
            logger(log, f'\t---> Genome coverage: {coverage}\n', mode="simple")
        
        # Assemble the genomes of the samples that passed the quality control
        assembly_dir = os.path.join(outDir, "assemblies")
        if qc_verdict == "Pass" and cov_verdict == "Pass":
            process_data.assemble(sample=sample, reads=fastqs, assembly_dir=assembly_dir, assembler=args.assembler, sequencer="Illumina", cpus=cpus, logfile=log)

            genome = os.path.join(assembly_dir, "genomes", f"{sample}.fasta")

            try:
                hit, tax_confirm, possibilities = find_organism.find_species_with_kmrf(s_name=sample, lab_species=sys_organism, genome=genome, dataOut=outDir, org_type="bacteria", logfile=log)
            except Exception as e:
                time_print(f"Error at find organism step: {e}", "Fail")
                logger(log, f"Error at find organism step: {e}")
                continue

            #######
            best_org = None
            best_percent = None
            best_other_org = None
            if tax_confirm == "Pass" or tax_confirm == "N/A":
                best_org = hit.split(" (")[0]
                best_percent = f"{float(hit.split(' (')[1].strip(')%')):.2f}"
            else:
                best_other_hit = hit.split(" --- ")[0].split(": ")[1]
                best_other_org = best_other_hit.split(" (")[0]
                best_other_percent = best_other_hit.split(" (")[1].strip(")")
            # Update coverage for samples with pre-unknown organisms
            if coverage == "N/A" and organism != "unknown":
                if best_org:
                    identified_org = best_org
                elif best_other_org:
                    identified_org = best_other_org
                if identified_org in bacteria:
                    g_size = bacteria.get(identified_org, None)[0]
                    if os.path.exists(qc_file):
                        with open(qc_file, 'r') as qcFile:
                            lines = qcFile.readlines()

                        bases_line = next((l for l in lines if "Total bases:" in l), None)
                        if bases_line is None:
                            bases_line = next((l for l in lines if "Tot bases:" in l), None)

                        total_bases = int(float(bases_line.split(":")[1].strip().replace(",", "")))

                        coverage = total_bases / int(g_size)
                        cov_display = f"{coverage:.2f}X"

            if tax_confirm == "Fail":
                writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, f'Closest: {best_other_org}', best_other_percent, f'{cov_display}', min_cov, cov_verdict, tax_confirm])
            else:
                writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, min_cov, cov_verdict, tax_confirm])

        else:
            best_org = "N/A"
            best_percent = "N/A"
            tax_confirm = "N/A"

            writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, min_cov, cov_verdict, tax_confirm])
# Assess Genome Quality with CheckM
genomes_dir = os.path.join(outDir, "assemblies", "genomes")
checkm_dir = os.path.join(outDir, "checkM")
checkm_out = process_data.checkM_stats(genomes_dir=genomes_dir, cpus=cpus, outdir=checkm_dir, logfile=log)

# Update the summary file with CheckM results
# qc_summary = os.path.join(qc_out, "qc_summary.tsv")
qc_summary = os.path.join(outDir, f"{run_name}_qc_summary.tsv")
time_print('Final summary', "Header")
logger(log, 'Final summary', "Header")

process_data.make_summary(qc_summary=qc_summary, temp_qc_summary=temp_qc_summary, header=header, checkm_out=checkm_out, logfile=log)

# Remove temp_dir
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
    time_print(f"Temporary directory removed.")
    logger(log, f"Temporary directory removed.")

# End the timer
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the execution time
end_message = f"PIPELINE COMPLETED: Processed {sample_number} samples in {int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"
time_print(end_message, "Header")
logger(log, end_message, "Header")
