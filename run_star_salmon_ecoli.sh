#!/bin/bash
#SBATCH --job-name=star_salmon
#SBATCH --mail-user=cfvall@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --account= FILAAAAAAAAAAAAAAAAAAAAAAA
#SBATCH --partition=standard
#SBATCH --output=/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/taxonomic_classification/star_salmon_ecoli/%x-%j_%a.log
#SBATCH --array=0-10  # Will be updated dynamically

# Exit on error
set -e

# Function for logging messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check if a directory exists
check_directory() {
    if [ ! -d "$1" ]; then
        log_message "ERROR: Directory $1 does not exist!"
        exit 1
    fi
}

# Function to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        log_message "ERROR: File $1 does not exist!"
        exit 1
    fi
}

# Load conda environment
log_message "Activating conda environment"
eval "$(conda shell.bash hook)"
conda activate general_bio || {
    log_message "ERROR: Failed to activate conda environment"
    exit 1
}

# Define paths
STAR_INDEX="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GCA_003697165_ecoli/STAR_index"
ECOLI_TRANSCRIPTS="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GCA_003697165_ecoli/GCF_000005845.2_ASM584v2_rna_from_genomic.fna"
READS_DIR="/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/fastq_files/11139/JB_1"
OUTPUT_DIR="/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/taxonomic_classification/star_salmon_ecoli"
THREADS=$SLURM_CPUS_PER_TASK

# Verify required directories and files exist
log_message "Checking input directories and files"
check_directory "$STAR_INDEX"
check_file "$ECOLI_TRANSCRIPTS"
check_directory "$READS_DIR"
check_directory "$OUTPUT_DIR"

# Generate sample list dynamically
SAMPLES=($(ls ${READS_DIR}/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
TOTAL_SAMPLES=${#SAMPLES[@]}

# Validate SLURM array range
if [ "$SLURM_ARRAY_TASK_ID" -ge "$TOTAL_SAMPLES" ]; then
    log_message "ERROR: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID exceeds available samples ($TOTAL_SAMPLES)"
    exit 1
fi

# Select sample for this job
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
STAR_OUT="${OUTPUT_DIR}/${SAMPLE_NAME}_star"
SALMON_OUT="${OUTPUT_DIR}/${SAMPLE_NAME}_salmon_quant"

# Create output directories
log_message "Creating output directories for sample $SAMPLE_NAME"
mkdir -p "$STAR_OUT"
mkdir -p "$SALMON_OUT"

# Check input files exist
R1="${READS_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz"
R2="${READS_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"
check_file "$R1"
check_file "$R2"

# Run STAR alignment
log_message "Starting STAR alignment for sample $SAMPLE_NAME"
STAR --runThreadN $THREADS \
     --genomeDir "$STAR_INDEX" \
     --readFilesIn "$R1" "$R2" \
     --readFilesCommand zcat \
     --alignIntronMax 1 \
     --outFileNamePrefix "${STAR_OUT}/" \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 6000000000 \
     2>&1 | tee "${STAR_OUT}/star.log"

# Check STAR output
if [ ! -f "${STAR_OUT}/Aligned.sortedByCoord.out.bam" ]; then
    log_message "ERROR: STAR alignment failed to produce BAM output"
    exit 1
fi

# Run Salmon quantification in alignment-based mode
log_message "Starting Salmon quantification for sample $SAMPLE_NAME"
salmon quant -t "$ECOLI_TRANSCRIPTS" -l A \
             -a "${STAR_OUT}/Aligned.sortedByCoord.out.bam" \
             -p $THREADS \
             --validateMappings \
             --gcBias \
             --seqBias \
             -o "$SALMON_OUT" \
             2>&1 | tee "${SALMON_OUT}/salmon.log"

# Check Salmon output
if [ ! -f "${SALMON_OUT}/quant.sf" ]; then
    log_message "ERROR: Salmon quantification failed to produce output file"
    exit 1
fi

log_message "Processing completed successfully for sample $SAMPLE_NAME"
