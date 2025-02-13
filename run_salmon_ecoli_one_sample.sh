#!/bin/bash
#SBATCH --job-name=salmon_ecoli
#SBATCH --mail-user=cfvall@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m  # 4GB per CPU (Total: 16GB)
#SBATCH --time=02:00:00  # Increased time limit to 2 hours
#SBATCH --account=prensner0
#SBATCH --partition=standard
#SBATCH --output=/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/taxonomic_classification/salmon_map_ecoli/%x-%A_%a.log
#SBATCH --array=0-5  # Adjust dynamically based on the number of samples

# Activate the correct environment
conda activate general_bio  # Ensure Salmon is installed

# Define paths
SALMON_INDEX="/nfs/turbo/umms-prensnerturbo/shared/assemblies/GCA_003697165_ecoli/salmon_ecoli_index"
READS_DIR="/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/fastq_files/11139"
OUTPUT_DIR="/scratch/prensner_root/prensner0/cfvall/CUTandRUN_Venneti_lab_2025/taxonomic_classification/salmon_map_ecoli"
THREADS=4

# Generate a list of sample names dynamically
SAMPLES=($(ls ${READS_DIR}/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))

# Select the sample based on the SLURM_ARRAY_TASK_ID
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Run Salmon quantification
salmon quant -i $SALMON_INDEX \
             -l A \
             -1 ${READS_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz \
             -2 ${READS_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz \
             -p $THREADS \
             --validateMappings \
             -o ${OUTPUT_DIR}/${SAMPLE_NAME}_salmon_quant

echo "Salmon alignment to E. coli completed for $SAMPLE_NAME."
