#!/bin/bash
#SBATCH --job-name=pipelineWrapper     # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email_address>    # Where to send mail	
#SBATCH --ntasks=16                    # Run on 16 cpus
#SBATCH --nodes=1                      # Run on a single node (computer)
#SBATCH --output=serial_test_%j.out    # Standard output and error log
#SBATCH --mem=32G                      # Allocated 32 gigabytes of memory for the job.

# Run this command from inside of the scripts directory. The format should look like this:
#<main folder>/
#├─ arriberelab_docker/
#│  ├─ *arribere lab scripts*/
#│  ├─ *THIS SLURM SCRIPT*
#│  ├─ RNASeqSettings.txt (edited for your run!)
#├─ data_dir/
#│  ├─ <unziped fastq files>.fastq
#├─ genome_dir/
#│  ├─ <mouse genome - indexed with STAR>
#├─ output_dir/
#│  ├─ <empty until after run>


# First set the locations of your genome and data and output files:
GENOME="../genome_dir"
DATA="../data_dir"  # Make sure that all the fastq files are already unzipped in here
OUTPUT="../output_dir"

# Load modules on the hummingbird
module load python/3.7
module load star/star2-2.5.3a
module load samtools
module load bcftools

# Iterate through all of the different fastq files, running the pipeline for each:
for fastq in $DATA/*.fastq ; do
  echo $fastq
  name="$(basename $fastq)"
  name="${name:0:-6}"
  echo ; echo "Starting to run pipeline for $name at $(date +"%D %H:%M")"
  mkdir $OUTPUT/$name
  python3 pipelineWrapperRNASeq.py $fastq RNASeqSettings.txt $OUTPUT"/"$name"/"$name
done
