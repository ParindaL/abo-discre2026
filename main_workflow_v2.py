import subprocess
import os

# ACTIVATE ENV WITH BIOINFO TOOLS 
# --- Configuration ---
# TODO: MODIFY THESE VARIABLES FOR YOUR ANALYSIS, CHECK FILE DIRECTORY
# -----------------------------------------------
print('CAUTION: BIOINFORMATICS ENVIRONMENT REQUIRED')
FASTQ_FILE = input("\nENTER FASTQ FILE NAME (test.fastq) :")
REFERENCE_GENOME = "hg38.fa"
OUTPUT_DIR = input("\nENTER OUTPUT FOLDER (automate_test) :") # output folder - will be created if does not exist.
BED_FILE = "designed.hg38.bed"  

# --- Helper Function to Run Commands ---
def run_command(command):
    """Prints and runs a command, checking for errors."""
    print(f"--- Running command --- \n{command}\n")
    try:
        subprocess.run(command, shell=True, check=True, text=True)
        print("--- Command successful ---\n")
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR: Command failed with exit code {e.returncode} ---")
        print(f"--- Failing command: {command} ---")
        # Exit the script if a command fails
        exit(1)

# --- Main Pipeline Workflow ---
def main():
    # --- 0. Setup ---
    print("--- Starting the bioinformatics pipeline ---")
    
    # Create the output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Derive a base name for all output files from the FASTQ
    sample_name = os.path.basename(FASTQ_FILE).split('.')[0]
    
    # Define file paths that will be created during the pipeline
    trimmed_fastq = os.path.join(OUTPUT_DIR, f"{sample_name}_trimmed.fq")
    sam_file = os.path.join(OUTPUT_DIR, f"{sample_name}.sam")
    bam_file = os.path.join(OUTPUT_DIR, f"{sample_name}.bam")
    sorted_bam = os.path.join(OUTPUT_DIR, f"{sample_name}.sorted.bam")
    final_bam = os.path.join(OUTPUT_DIR, f"{sample_name}.final.bam")
    raw_vcf = os.path.join(OUTPUT_DIR, f"{sample_name}.raw.vcf")
    annotated_vcf = os.path.join(OUTPUT_DIR, f"{sample_name}.annotated.vcf")

    # --- 1. Quality Control & Trimming ---
    print("--- Step 1: Running FastQC and Trim Galore! ---")

    run_command(f"fastqc {FASTQ_FILE} -o {OUTPUT_DIR}")
    run_command(f"trim_galore --output_dir {OUTPUT_DIR} {FASTQ_FILE}")

    # --- 2. Alignment ---
    print("--- Step 2: Aligning reads with BWA-MEM ---")

    base_trimmed_fastq = os.path.join(OUTPUT_DIR, f"{sample_name}_trimmed.fq")
    run_command(f"bwa mem {REFERENCE_GENOME} {base_trimmed_fastq} > {sam_file}")
    
    # --- 3. Post-Alignment Processing ---
    print("--- Step 3: Processing SAM/BAM file ---")
    run_command(f"samtools view -bS {sam_file} > {bam_file}")
    run_command(f"samtools sort {bam_file} -o {sorted_bam}")

    # ---  Add Read Group Information ---
    print("--- Step 3b: Adding Read Group information ---")
    rg_bam = os.path.join(OUTPUT_DIR, f"{sample_name}.rg.bam")
    run_command(
        f"gatk AddOrReplaceReadGroups "
        f"-I {sorted_bam} "
        f"-O {rg_bam} "
        f"-RGID 4 "
        f"-RGLB lib1 "
        f"-RGPL IONTORRENT "
        f"-RGPU unit1 "
        f"-RGSM {sample_name}"  # Use the sample name derived from the FASTQ
    )
    
    # --- Now use the NEW rg_bam file for MarkDuplicates ---
    # run_command(f"gatk MarkDuplicates -I {rg_bam} -O {final_bam} -M {OUTPUT_DIR}/dup_metrics.txt")
    run_command(f"samtools index {rg_bam}")
    
    # --- 4. Variant Calling ---
    print("--- Step 4: Calling variants with GATK HaplotypeCaller ---")
    run_command(f"gatk HaplotypeCaller -R {REFERENCE_GENOME} -I {rg_bam} -O {raw_vcf} ") #-L {BED_FILE}
   
    print("--- Pipeline finished successfully! ---")
    print("--- You can now proceed to variant annotation. ---")

# --- Run the main function when the script is executed ---
if __name__ == "__main__":
    main()
