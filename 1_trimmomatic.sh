#!/bin/bash

#BATCH --job-name=trimmomatic #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=10G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


##########################
# Module load:
module load java/x64/17u3

#########################
# Run commands from here:

input_dir="/groups/envgenom/Thais_DalSasso/RNAseq/resequencing_Zpa/rawdata"
output_dir="/home/dalsasso/RNAseq/trimmed/resequencing_Zpa/trimmed"


trimmomatic_jar="/data/biosoftware/Trimmomatic/Trimmomatic/trimmomatic-0.39.jar"
adapters_file="/home/dalsasso/RNAseq/rawdata/adapters.fa"

# Number of threads
threads=15


mkdir -p "$output_dir"

# Loop through all R1 files in the input directiry
for infile in "$input_dir"/*_R1_001.fastq.gz; 
do
    base=$(basename "${infile}" _R1_001.fastq.gz)
    infile_r2="${base}_R2_001.fastq.gz"

    # Run Trimmomatic
    java -jar "$trimmomatic_jar" PE -threads "$threads" \
        "${infile}" "${input_dir}/${infile_r2}" \
        "${output_dir}/${base}_R1.paired.trimmed.fastq.gz" "${output_dir}/${base}_R1.unpaired.trimmed.fastq.gz" \
        "${output_dir}/${base}_R2.paired.trimmed.fastq.gz" "${output_dir}/${base}_R2.unpaired.trimmed.fastq.gz" \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:"$adapters_file":2:40:15 

done
