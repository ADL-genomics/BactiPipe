import os
import subprocess
from bactipipe.scripts.utils import time_print, simple_print, logger
from .config import DATA_DIR

kmerfinder_version = "3.0.2" # To be updated manually

def parse_kkn(kkn_report, logfile=None):
    log = logfile
    time_print("Parsing Kraken2 report...")
    logger(log, "Parsing Kraken2 report...")
    with open(kkn_report, 'r') as ts:
        species = {}
        for line in ts:
            perc_id, frag1, frag2, level, tax_id, name = line.strip().split("\t")
            if level == "S":
                species[tax_id] = [float(perc_id.strip()), name.strip()]
                if float(perc_id) >= 0.1:
                    print(f'\t----> {tax_id}\t{perc_id.strip()}\t{name.strip()}')
        return species

def find_species_with_kkn(reads, kkn_db, s_name, expected_taxonomy, dataOut, org_type="bacteria", logfile=None):
    log = logfile
    intro = f"Running Kraken2 analysis for species confirmation ..."
    time_print(intro)
    logger(log, intro)

    out_folder = os.path.join(dataOut, f'kraken_out/{s_name}')
    lab_species = expected_taxonomy[1]
    lab_taxid = expected_taxonomy[0]

    # Output folder
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    kkn_report = f"{out_folder}/{s_name}_kkn.report"


    # Running the kraken2 command
    vercmd = f'kraken2 --version'
    stdout = subprocess.run(vercmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = stdout.stdout.decode('utf-8')
    for line in out.split('\n'):
        if line.startswith('Kraken 2'):
            version = line.split()[-1]
            break

    kkninfo = f"Running Kraken2 (version {version}) analysis for species confirmation on {s_name}..."
    time_print(kkninfo)
    logger(log, kkninfo)

    if len(reads) == 2:
        cmd = f"kraken2 --db {kkn_db} --output {out_folder}/{s_name}_kkn.out --report {kkn_report} --paired {reads[0]} {reads[1]}"

    elif len(reads) == 1:
        cmd = f"kraken2 --db {kkn_db}  --output {out_folder}/{s_name}_kkn.out --report {kkn_report} {reads[0]}"
    
    if not os.path.exists(kkn_report):
        cmdinfo = f"Running: {cmd}"
        time_print(cmdinfo, s_type='command')
        logger(log, cmdinfo, s_type='command')
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
    else:
        info = f"Skipping Kraken2 analysis. Report already exists"  
        time_print(info)
        logger(log, info)

    # Parsing the results
    specs = []
    hits = []
    other_displ = []
    best_hit = "N/A"
    best_taxid = "N/A"
    best_perc_id = 0
    max_perc_id = 0
    spec_found = False

    kkn_out = parse_kkn(kkn_report)
    
    for taxid in kkn_out.keys():
        perc_id, name = kkn_out[taxid]
        if perc_id >= 10: # Minimum percentage identity (will be a parameter)
            spec_found = True
            if perc_id > best_perc_id:
                if best_taxid != "N/A":
                    other_displ.append(f'{best_taxid}:{best_hit} ({best_perc_id})')
                best_taxid = taxid
                best_perc_id = perc_id
                best_hit = name
            else:
                other_displ.append(f'{taxid}:{name} ({perc_id})')
        else:
            if perc_id > max_perc_id:
                max_perc_id = perc_id
    
            other_displ.append(f'{taxid}:{name} ({perc_id})')

    if spec_found == True:
        bh_display = f'{best_taxid}:{best_hit} ({best_perc_id})'
    else:
        bh_display = f'{best_taxid}:{best_hit} ({max_perc_id})'

    if len(other_displ) > 0:
        others = ', '.join(other_displ)
    else:
        others = 'N/A'


    if best_taxid == lab_taxid:
        tax_confirm = 'Pass'
    else:
        tax_confirm = 'Fail'
        bh_display = f'Determined: {bh_display}) --- Expected: {lab_taxid}:{lab_species}'

    return bh_display, tax_confirm, others


def find_species_with_kmrf(s_name, lab_species, genome, dataOut, org_type="bacteria", logfile=None):
    log = logfile
    intro = f"Sample {s_name}: Running KmerFinder (version: {kmerfinder_version})"
    out_folder = os.path.join(dataOut, f'kmerFinder/{s_name}')

    # Output folder
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Database files
    db_dir = DATA_DIR
    kmer_db = os.path.join(db_dir, "kmerfinder_db", org_type, f"{org_type}.ATG")
    taxa_file = os.path.join(db_dir, "kmerfinder_db", org_type, f"{org_type}.tax")

    # Running the kmerfinder command
    cmd = f"kmerfinder.py -i {genome} -o {out_folder} -db {kmer_db} -tax {taxa_file}"
    results_file = os.path.join(out_folder, "results.txt")

    if not os.path.exists(results_file):
        time_print(intro)
        logger(log, intro)

        cmdinfo = f"--> Command: {cmd}"
        simple_print(cmdinfo, s_type='command')
        logger(log, cmdinfo, s_type='command', mode="simple")
        try:
            process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            stdout = out.decode('utf-8', errors='replace')
            stderr = err.decode('utf-8', errors='replace')
            # logger(log, f"KmerFinder stdout:\n{stdout}")
            if stderr:
                logger(log, f"KmerFinder stderr:\n{stderr}")
            if process.returncode != 0:
                err_info = f"KmerFinder failed for sample {s_name}: exit {process.returncode}"
                time_print(err_info, "Fail")
                logger(log, err_info)
                return
            else:
                success_info = "Command exit status: Success!"
                time_print(success_info, "Pass")
                logger(log, success_info)
        except Exception as e:
            err_msg = f"Error running KmerFinder for sample {s_name}: {e}"
            time_print(err_msg, "Fail")
            logger(log, err_msg)
            return
    else:
        info = f"Sample {s_name}: Skipping KmerFinder. Results already exist."
        time_print(info)
        logger(log, info)

    # Parsing the results
    unsorted_hits = {}
    other_displ = []
    with open(results_file, "r") as file:
        for row in file:
            if not row.startswith("#"):
                parts = row.strip().split('\t')
                description = parts[14].strip()
                species = parts[-1]
                if species == 'Streptococcus equi' and 'zooepidemicus' in description:
                    species = 'Streptococcus equi subsp. zooepidemicus'
                elif species == 'Streptococcus equi':
                    species = 'Streptococcus equi subsp. equi'
                qcov = round(float(parts[5].strip()), 2)
                unsorted_hits[species] = unsorted_hits.get(species, 0) + qcov

    hits = sorted(unsorted_hits.items(), key=lambda item: item[1], reverse=True)
    best_hit = hits[0][0]
    bh_display = f'{hits[0][0]} ({hits[0][1]})'

    if len(hits) > 1:
        other_displ = [f"{hit[0]} ({hit[1]:.2f})" for hit in hits[1:]]
        others = ', '.join(other_displ)
    else:
        others = 'N/A'
    
    if lab_species.lower() == "unknown":
        tax_confirm = 'N/A'
    elif best_hit == lab_species:
        tax_confirm = 'Pass'
    else:
        tax_confirm = 'Fail'
        bh_display = f"Determined: {bh_display}) --- Expected: {lab_species}"

    return bh_display, tax_confirm, others
