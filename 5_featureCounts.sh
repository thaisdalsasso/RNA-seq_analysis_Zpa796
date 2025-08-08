#!/bin/bash

#BATCH --job-name=featureCounts #Give your job a name.

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
#module load $MODULE_NAME 

#########################
# Run commands from here:


#Set input files and strings
work_dir="/home/dalsasso/RNAseq/4_mapping2fungus/6448_resequencing_Zpa"
alignment_dir="/home/dalsasso/RNAseq/4_mapping2fungus/6448_resequencing_Zpa/alignments"
annotation_file="/home/dalsasso/data/References/gene_models/Zpa796.braker.rna.bbd2.gff3"
create_dir="counts"
annotation_type="mRNA" #i.e: CDS, gene, exon (default)...
gff_atribute="Parent" #i.e.: gene_id (default), transcript_id...


if [ ! -d "$work_dir/$create_dir" ]; then
    mkdir -p "$work_dir/$create_dir"
fi

cd "$work_dir/$create_dir"


/home/dalsasso/softwares/subread-2.0.6-Linux-x86_64/bin/featureCounts -s 0 -p -T 15 \
-t "$annotation_type" -g "$gff_atribute" -a "$annotation_file" \
-o ./Quantified_Zpa796_project6448_all_invitro_inplanta.txt \
$(ls "$alignment_dir/mapped2Zpa_invitro"/Zpa796_*bam) \
$(ls "$alignment_dir/mapped2Zpa_inplanta"/Zpa796_HmZpa796*bam)














