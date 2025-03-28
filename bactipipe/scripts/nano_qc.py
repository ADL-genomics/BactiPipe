import os
import subprocess
from bactipipe.scripts.utils import time_print

def qc_nano(
    fastq_file=None,
    raw_folder=None,
    genome_size=5e6,
    min_avg_quality=15,
    min_coverage=100,
    desired_coverage=300,
    step_coverage=50,
    output_fastq="trimmed_reads.fastq.gz",
    output_dir=".",
    cpus="8",
    single=True
):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_fastq = os.path.join(output_dir, output_fastq)

    if raw_folder:
        fastq_file = os.path.join(output_dir, os.path.basename(raw_folder) + ".fastq.gz")
        if single:
            time_print(f"--->Processing: {os.path.basename(raw_folder)}--->{os.path.basename(output_fastq)}\n")
            time_print("Consolidating nanopore raw reads into one file...")
        if not os.path.exists(fastq_file):
            subprocess.run(f"cat {raw_folder}/*.fastq.gz > {fastq_file}", shell=True)
        else:
            if single:
                time_print(f"File {fastq_file} already exists. Skipping consolidation.")
    else:
        if single:
            time_print(f"--->Processing: {os.path.basename(fastq_file)}--->{os.path.basename(output_fastq)}\n")

    # Step 1: Run NanoPlot on the raw fastq file to obtain quality metrics
    final_qual_file = f"{output_dir}/quality_metrics_after_qc.txt"
    if not os.path.exists(final_qual_file):
        if single:
            time_print("Running NanoPlot to assess quality metrics...")
        initial_out = os.path.join(output_dir, "nanoplot_initial")
        nanoplot_command = [
            "NanoPlot", "--fastq", fastq_file,
            "--outdir", initial_out,
            "--threads", cpus, "--plots", "dot"
        ]
        subprocess.run(nanoplot_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Parse the NanoPlot summary file
        nanoplot_summary_file = os.path.join(initial_out, "NanoStats.txt")
        with open(nanoplot_summary_file, 'r') as f:
            lines = f.readlines()

        # Extract relevant data
        total_bases = int(float(next(line for line in lines if "Total bases:" in line).split(":")[1].strip().replace(",", "")))
        if single:
            time_print(f"\t---> Total bases: {total_bases:,}")
        mean_quality = float(next(line for line in lines if "Mean read quality:" in line).split(":")[1].strip())
        if single:
            time_print(f"\t---> Mean quality: {mean_quality}")

        # Check if there is enough data for at least 100X coverage
        min_required_bases = min_coverage * genome_size
        if total_bases < min_required_bases:
            if single:
                print(f"Insufficient data for {min_coverage}X coverage. Total bases available: {total_bases:,}. Required: {min_required_bases:,}.")
            os.rename(os.path.join(initial_out, "LengthvsQualityScatterPlot_dot.png"), f"{output_dir}/raw_reads_quality.png")
            os.rename(os.path.join(initial_out,"NanoStats.txt"), f"{output_dir}/raw_reads_quality_metrics.txt")
            return

        # Step 2: Iteratively trim the reads
        coverage_target = desired_coverage
        while coverage_target >= min_coverage:
            # Calculate the required number of bases for the target coverage
            required_bases = int(coverage_target * genome_size)

            # Check if there is enough data for the current target coverage
            if total_bases < required_bases:
                if single:
                    time_print(f"Insufficient data for {coverage_target}X coverage. Total bases available: {total_bases}, Required: {required_bases:,}.")
                coverage_target -= step_coverage
                continue

            # Trim reads using filtlong
            if single:
                time_print(f"\nFiltering reads to achieve {coverage_target}X coverage...")
            filtlong_command = ["filtlong", "--target_bases", str(required_bases), fastq_file]
            trimmed_output = os.path.join(output_dir, f"trimmed_{coverage_target}X.fastq")
            with open(trimmed_output, "w") as trimmed_file:
                subprocess.run(filtlong_command, stdout=trimmed_file, stderr=subprocess.DEVNULL)
            if single:
                time_print("Compressing trimmed reads...")
            gzip_command = ["gzip", "-f", trimmed_output]
            subprocess.run(gzip_command, check=True)

            # Recalculate quality metrics on the trimmed reads
            if single:
                print(f"Running NanoPlot on trimmed reads for {coverage_target}X coverage...")
            last_out = os.path.join(output_dir, "nanoplot_last")
            nanoplot_command = [
                "NanoPlot", "--fastq", f"{trimmed_output}.gz",
                "--threads",  cpus, "--plots", "dot",
                "--outdir", last_out
            ]

            subprocess.run(nanoplot_command)

            # Parse new quality metrics
            summary_last = os.path.join(last_out, "NanoStats.txt")
        
            with open(summary_last, "r") as np_summary:
                np_lines = np_summary.readlines()
            
            new_mean_quality = float(next(line for line in np_lines if "Mean read quality:" in line).split(":")[1].strip())

            # Check if the new quality meets the threshold
            if new_mean_quality >= min_avg_quality:
                if single:
                    time_print("\nQuality threshold met. Saving data...\n")
                os.rename(f"{trimmed_output}.gz", output_fastq)
                os.rename(os.path.join(initial_out, "LengthvsQualityScatterPlot_dot.png"), f"{output_dir}/raw_reads_quality.png")
                os.rename(os.path.join(last_out,"LengthvsQualityScatterPlot_dot.png"), f"{output_dir}/quality_after_qc.png")
                os.rename(os.path.join(last_out, "NanoStats.txt"), final_qual_file)
                os.rename(os.path.join(initial_out,"NanoStats.txt"), f"{output_dir}/raw_reads_quality_metrics.txt")

                if single:
                    time_print(f"Trimmed reads saved to {output_fastq} with {coverage_target}X coverage and average quality {new_mean_quality:.2f}.")
                return
            else:
                os.rename(os.path.join(initial_out, "LengthvsQualityScatterPlot_dot.png"), f"{output_dir}/raw_reads_quality.png")
                os.rename(os.path.join(initial_out,"NanoStats.txt"), f"{output_dir}/raw_reads_quality_metrics.txt")
                if single:
                    time_print(f"Quality {new_mean_quality:.2f} below threshold for {coverage_target}X coverage.")
                os.remove(f"{trimmed_output}.gz")
                coverage_target -= step_coverage
        
        # If we exit the loop without meeting the requirements
        if single:
            time_print(f"Unable to achieve {min_avg_quality} average quality with minimum coverage of {min_coverage}X.")

        # Cleanup
        if single:
            time_print("Cleaning up intermediate files...")
        if os.path.exists(initial_out):
            subprocess.run(["rm", "-r", initial_out])
        if os.path.exists(last_out):
            subprocess.run(["rm", "-r", last_out])

        if coverage_target < min_coverage:
            print(f"Unable to achieve {min_avg_quality} average quality with minimum coverage of {min_coverage}X.")

    else:
        if single:
            time_print("Quality metrics file already exists. Skipping quality assessment.")
