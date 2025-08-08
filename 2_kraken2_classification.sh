#!/bin/bash

#BATCH --job-name=kraken2 #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=20G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


##########################
# Module load:
#module load $MODULE_NAME 

#########################
# Run commands from here:


work_dir="/home/dalsasso/RNAseq/trimmed/resequencing_Zpa/trimmed"
database_dir="/data/biosoftware/kraken2/db/minikraken_8GB_20200312"

cd "$work_dir"

mkdir -p kraken2
mkdir -p kraken2/classified
mkdir -p kraken2/unclassified

for input in *_R1.paired.trimmed.fastq.gz; do
    input2=$(echo "$input" | sed "s/R1.paired/R2.paired/g")
    output1=$(echo "$input" | sed "s/_R1.paired.trimmed.fastq.gz//g")

    kraken2 --use-names --threads 20 --gzip-compressed \
    --db "$database_dir" \
    --report ./kraken2/"$output1"_report.txt \
    --paired ./"$input" ./"$input2" \
    --classified-out ./kraken2/classified/"$output1"_paired.trimmed.contamination#.fastq \
    --unclassified-out ./kraken2/unclassified/"$output1"_paired.trimmed.filtered#.fastq \
    > ./kraken2/"$output1"_kraken.txt

    gzip ./kraken2/classified/"$output1"_paired.trimmed.contamination*.fastq  \
    gzip ./kraken2/unclassified/"$output1"_paired.trimmed.filtered*.fastq \

done








