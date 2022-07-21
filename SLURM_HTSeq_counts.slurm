#!/bin/bash
#SBATCH --job-name=HTSeq-counts        # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mviscard@ucsc.edu    # Where to send mail	
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1                      # Run on a single node (computer)
#SBATCH --output=serial_test_%j.out    # Standard output and error log
#SBATCH --mem=32G                      # Allocated 32 gigabytes of memory for the job.

pip3 install htseq

python3 -h HTSeq.scripts.count /hb/home/saguiar/220721_RNASeq/output_dir/*/*.bam /hb/home/saguiar/220721_RNASeq/genome_dir/*.gtf
