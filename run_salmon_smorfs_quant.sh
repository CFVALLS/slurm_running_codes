#!/bin/bash

# Set strict error handling
set -euo pipefail
IFS=$'\n\t'

# Script configuration
readonly BAM_LIST="bam_list.txt"
readonly OUTPUT_DIR="salmon_output"
readonly GTF_FILE="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/gencode.v45.Riboseq-ORF.with_chr.gtf"
readonly GENOME_FA="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GENCODE_45/GRCh38.primary_assembly.genome.fa"
readonly LOG_DIR="logs"
readonly MAX_ATTEMPTS=2

# SLURM settings
readonly CPUS=8
readonly MEMORY=16G
readonly TIME="12:00:00"
readonly EMAIL="cfvall@umich.edu"
readonly ACCOUNT=""  # This should be set to your SLURM account

# Function to log messages with timestamps
log() {
    local message="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" | tee -a "${LOG_DIR}/pipeline.log"
}

# Function to check if required files exist
check_requirements() {
    local missing_files=0

    for file in "${BAM_LIST}" "${GTF_FILE}" "${GENOME_FA}"; do
        if [[ ! -f "${file}" ]]; then
            log "Error: Required file '${file}' not found."
            ((missing_files++))
        fi
    done

    # Check if Salmon is available
    if ! command -v salmon &> /dev/null; then
        log "Error: Salmon command not found. Please ensure it's installed and in PATH."
        ((missing_files++))
    fi

    # Check if gffread is available
    if ! command -v gffread &> /dev/null; then
        log "Error: gffread command not found. Please ensure it's installed and in PATH."
        ((missing_files++))
    fi

    if ((missing_files > 0)); then
        exit 1
    fi
}

# Function to create required directories
setup_directories() {
    mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
    log "Created output directories: ${OUTPUT_DIR} and ${LOG_DIR}"
}

# Function to activate the required environment
activate_environment() {
    # Save current IFS and restore it after source
    local OLD_IFS="$IFS"
    IFS=$' \t\n'
    
    # Try multiple paths for conda activation
    if source activate general_bio 2>/dev/null || \
       source ~/miniconda3/bin/activate general_bio 2>/dev/null || \
       source /usr/local/miniconda3/bin/activate general_bio 2>/dev/null; then
        log "Successfully activated conda environment"
    else
        log "Failed to activate conda environment"
        exit 1
    fi
    
    # Restore IFS
    IFS="$OLD_IFS"
}

# First slurm job to run gffread
submit_gffread_job() {
    local cds_file="$1"
    local cds_dir
    cds_dir="$(dirname "${cds_file}")"

    # Create directory if it doesn't exist
    mkdir -p "${cds_dir}"

    # Check if the CDS file already exists and is non-empty
    if [[ -f ${cds_file} && -s ${cds_file} ]]; then
        log "CDS file '${cds_file}' already exists and is non-empty. Skipping extraction."
        return 0
    fi

    log "Extracting CDS sequences from '${GENOME_FA}'"
    if ! gffread -x "${cds_file}" -g "${GENOME_FA}" "${GTF_FILE}" -E; then
        log "Error: gffread failed to extract CDS sequences."
        rm -f "${cds_file}"  # Remove potentially partial output
        exit 1
    fi

    # Verify the output file exists and is non-empty
    if [[ ! -s "${cds_file}" ]]; then
        log "Error: Generated CDS file is empty or was not created."
        exit 1
    fi

    log "CDS sequences extracted to '${cds_file}'"
}

# Function to create and submit SLURM job for a sample
submit_job() {
    local bam_file="$1"
    local sample_name
    sample_name=$(basename "${bam_file}" .Aligned.out.bam)
    local output_dir="${OUTPUT_DIR}/${sample_name}"
    local cds_file="${OUTPUT_DIR}/cds/${sample_name}_cds.fa"
    local slurm_script="${LOG_DIR}/slurm_${sample_name}.sh"
    local log_file="${LOG_DIR}/${sample_name}.log"

    # Verify BAM file exists and is not empty
    if [[ ! -s "${bam_file}" ]]; then
        log "Error: BAM file '${bam_file}' is empty or does not exist."
        return 1
    }

    # Skip if output already exists and contains expected files
    if [[ -d "${output_dir}" && -f "${output_dir}/quant.sf" ]]; then
        log "Sample ${sample_name} has already been processed and contains output files. Skipping."
        return 0
    fi

    # Ensure the CDS file exists
    submit_gffread_job "${cds_file}"

    # Create SLURM script
    cat > "${slurm_script}" <<EOL
#!/bin/bash
#SBATCH --job-name=salmon_${sample_name}
#SBATCH --output=${LOG_DIR}/${sample_name}.out
#SBATCH --error=${LOG_DIR}/${sample_name}.err
#SBATCH --time=${TIME}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
#SBATCH --account=${ACCOUNT}

# Exit on error
set -euo pipefail

# Function to clean up incomplete results
cleanup() {
    if [[ \$? -ne 0 ]]; then
        echo "Error occurred, cleaning up incomplete results..."
        rm -rf "${output_dir}"
    fi
}
trap cleanup EXIT

# Load environment
$(declare -f activate_environment)
activate_environment

# Create output directory
mkdir -p "${output_dir}"

# Step 2: Run Salmon quantification
echo "[$(date)] Running Salmon quantification for ${sample_name}"
salmon quant \
    -t "${cds_file}" \
    -l U \
    -a "${bam_file}" \
    -o "${output_dir}" \
    -p "${CPUS}" \
    --geneMap "${GTF_FILE}" \
    2>&1 | tee -a "${log_file}"

# Verify the output files exist
if [[ ! -f "${output_dir}/quant.sf" ]]; then
    echo "Error: Expected output file quant.sf not found."
    exit 1
fi

echo "[$(date)] Salmon quantification completed for ${sample_name}"
EOL

    # Submit job with retry logic
    local attempt=1
    local job_id=""
    
    while [[ ${attempt} -le ${MAX_ATTEMPTS} ]]; do
        log "Submitting job for ${sample_name} (Attempt ${attempt})"
        if job_id=$(sbatch "${slurm_script}" 2>&1 | grep -o '[0-9]*'); then
            log "Submitted job ${job_id} for sample ${sample_name}"
            echo "${job_id}:${sample_name}:${slurm_script}" >> "${LOG_DIR}/submitted_jobs.txt"
            return 0
        else
            log "Error: Failed to submit job for ${sample_name} (Attempt ${attempt})"
            sleep 5  # Wait before retrying
            ((attempt++))
        fi
    done

    log "Error: Reached maximum attempts for ${sample_name}. Job not submitted."
    return 1
}

main() {
    log "Starting Salmon processing pipeline"

    # Check requirements and setup
    check_requirements
    setup_directories

    # Create a separate directory for CDS files
    mkdir -p "${OUTPUT_DIR}/cds"

    # Process each BAM file
    local failed_submissions=0
    local total_samples=0
    
    while IFS= read -r bam_file || [[ -n "$bam_file" ]]; do
        ((total_samples++))
        if [[ -f "${bam_file}" ]]; then
            submit_job "${bam_file}" || ((failed_submissions++))
        else
            log "Warning: BAM file not found: ${bam_file}"
            ((failed_submissions++))
        fi
    done < "${BAM_LIST}"

    # Report completion status
    log "Pipeline completed. Processed ${total_samples} samples with ${failed_submissions} failures."
    if ((failed_submissions > 0)); then
        log "Warning: Some jobs failed to submit. Check logs in ${LOG_DIR}"
        exit 1
    fi
}

# Run main function
main
