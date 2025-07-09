#!/bin/bash

#SBATCH -J dry_run        # Job name
#SBATCH -p batch                  # Partition
#SBATCH -N 1                      # Number of nodes
#SBATCH --cpus-per-task=4         # cpu per task
#SBATCH --mem=8G                  # memory
#SBATCH -t 36:00:00               # Maximum runtime
#SBATCH -o logs/dry_run_%j.out  # Standard output log
#SBATCH -e logs/dry_run_%j.err  # Standard error log

source activate /home/groups/hoolock2/u0/bd/miniconda3/envs/snake

snakemake -n --use-conda --jobs 50 --latency-wait 60 --keep-going --rerun-incomplete --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr} -J {rule}.{wildcards}"
