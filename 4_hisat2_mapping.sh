#!/bin/bash

#BATCH --job-name=hisat2 #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15 #Multithreading.
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

work_dir="/home/dalsasso/RNAseq/4_mapping2fungus/6448_resequencing_Zpa/alignments"
genome_index="/home/dalsasso/data/References/genomes/Z796IPnp_renamed.fasta.index"
reads_dir="/home/dalsasso/RNAseq/3_rRNA_free/6448_resequencing_Zpa/inplanta"
create_dir="mapped2Zpa_inplanta"


if [ ! -d "$work_dir/$create_dir" ]; then
    mkdir -p "$work_dir/$create_dir"
fi

cd "$work_dir/$create_dir"

for input_file in "${reads_dir}"/*_paired.trimmed.filtered.rRNAfree_PE1.fastq.gz; do
    input_file2=$(echo "$input_file" | sed 's/rRNAfree_PE1/rRNAfree_PE2/')
    
    RG_id=$(basename "$input_file" | sed 's/_paired.*//')
    RG_name=$(basename "$input_file" | sed 's/_paired.*//')
    
    echo "Analyzing $input_file and $input_file2"

    hisat2 -p 12 --rg-id "$RG_id" --rg SM:"$RG_name" \
    --summary-file "${RG_name}_summary.txt" \
    -x "$genome_index" -1 "$input_file" -2 "$input_file2" \
    -S "./${RG_name}.sam"
    
    samtools sort -@ 15 "./${RG_name}.sam" > "./${RG_name}.bam"
    
    samtools index -@ 15 "./${RG_name}.bam"
    
    rm "./${RG_name}.sam"
done

