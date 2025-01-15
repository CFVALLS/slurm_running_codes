#!/bin/bash

# Set strict error handling
set -euo pipefail
IFS=$'\n\t'

# Script configuration
readonly BAM_LIST="bam_list.txt"
readonly OUTPUT_DIR="featurecounts_output"
readonly GTF_FILE="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/gencode.v45.Riboseq-ORF.with_chr.gtf"
readonly LOG_DIR="logs"
readonly MAX_ATTEMPTS=3

# SLURM settings
readonly CPUS=16
readonly MEMORY=4G
readonly TIME="2:00:00"
readonly EMAIL="cfvall@umich.edu"
readonly ACCOUNT=""

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check if required files exist
check_requirements() {
    if [[ ! -f "${BAM_LIST}" ]]; then
        log "Error: BAM list file '${BAM_LIST}' not found."
        exit 1
    fi

    if [[ ! -f "${GTF_FILE}" ]]; then
        log "Error: GTF file '${GTF_FILE}' not found."
        exit 1
    fi

    # Check if featureCounts is available
    if ! command -v featureCounts &> /dev/null; then
        log "Error: featureCounts command not found. Please ensure it's installed and in PATH."
        exit 1
    fi
}

# Function to create required directories
setup_directories() {
    mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
    log "Created output directories: ${OUTPUT_DIR} and ${LOG_DIR}"
}

# Function to create and submit SLURM job for a sample
submit_job() {
    local bam_file="$1"
    local sample_name
    sample_name=$(basename "${bam_file}" .Aligned.out.bam)
    local output_file="${OUTPUT_DIR}/${sample_name}_featurecounts.txt"
    local slurm_script="slurm_${sample_name}.sh"
    local log_file="${LOG_DIR}/${sample_name}.log"

    # Skip if output already exists
    if [[ -f "${output_file}" ]]; then
        log "Sample ${sample_name} has already been processed. Skipping."
        return 0
    fi

    # Create SLURM script
    cat > "${slurm_script}" <<EOL
#!/bin/bash
#SBATCH --job-name=fc_${sample_name}
#SBATCH --output=${LOG_DIR}/${sample_name}.out
#SBATCH --error=${LOG_DIR}/${sample_name}.err
#SBATCH --time=${TIME}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
#SBATCH --account=${ACCOUNT}

# Exit on error
set -e

# Load environment
source activate general_bio || source ~/miniconda3/bin/activate general_bio || {
    echo "Failed to activate conda environment"
    exit 1
}

# Run FeatureCounts with detailed logging
echo "[$(date)] Starting FeatureCounts for ${sample_name}"
featureCounts -T "${CPUS}" \
    -a "${GTF_FILE}" \
    -t CDS \
    -g gene_id \
    -o "${output_file}" \
    -p -B -C -M --fraction \
    "${bam_file}" 2>&1 | tee -a "${log_file}"

echo "[$(date)] Finished processing ${sample_name}"

# Verify output file was created
if [[ ! -f "${output_file}" ]]; then
    echo "Error: Output file was not created"
    exit 1
fi
EOL

    # Submit job with retry logic
    local attempt=1
    while [[ ${attempt} -le ${MAX_ATTEMPTS} ]]; do
        log "Submitting job for ${sample_name} (Attempt ${attempt})"
        local job_id
        job_id=$(sbatch "${slurm_script}" | grep -o '[0-9]*')

        if [[ -n "${job_id}" ]]; then
            log "Submitted job ${job_id} for sample ${sample_name}"
            echo "${job_id}" >> "${LOG_DIR}/submitted_jobs.txt"
            rm "${slurm_script}"  # Clean up SLURM script after successful submission
            return 0
        else
            log "Error: Failed to submit job for ${sample_name} (Attempt ${attempt})"
            ((attempt++))
        fi
    done

    log "Error: Reached maximum attempts for ${sample_name}. Job not submitted."
    return 1
}

main() {
    log "Starting FeatureCounts processing pipeline"
    
    # Check requirements and setup
    check_requirements
    setup_directories

    # Process each BAM file
    local failed_submissions=0
    while IFS= read -r bam_file; do
        if [[ -f "${bam_file}" ]]; then
            submit_job "${bam_file}" || ((failed_submissions++))
        else
            log "Warning: BAM file not found: ${bam_file}"
            ((failed_submissions++))
        fi
    done < "${BAM_LIST}"

    # Report completion status
    log "Pipeline completed. Failed submissions: ${failed_submissions}"
    if ((failed_submissions > 0)); then
        log "Warning: Some jobs failed to submit. Check logs in ${LOG_DIR}"
        exit 1
    fi
}

# Run main function
main
