#!/bin/bash
#SBATCH --job-name=HTSeq-counts        # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email_address>    # Where to send mail
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1                      # Run on a single node (computer)
#SBATCH --output=serial_test_%j.out    # Standard output and error log
#SBATCH --mem=32G                      # Allocated 32 gigabytes of memory for the job.

# First set the locations of your genome, and output directories:
GENOME="../genome_dir"
OUTPUT="../output_dir"
# Also make sure you change the email up above!!

pip3 install htseq

python3 -m HTSeq.scripts.count $OUTPUT/*/*.bam $GENOME/*.gtf
