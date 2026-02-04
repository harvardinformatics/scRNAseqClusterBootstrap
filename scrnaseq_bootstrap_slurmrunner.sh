#!/bin/bash
#SBATCH -p sapphire,shared
#SBATCH -e scrnaboot_%A.err
#SBATCH -o scrnaboot_%A.out
#SBATCH -J scrnaboot
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 8000
#SBATCH -t 36:00:00

module purge
module load python
conda activate my_snakemake

#which snakemake
#snakemake --version

snakemake --unlock --snakefile workflow/Snakefile --configfile config/config.yaml --use-conda --workflow-profile profiles/slurm --profile cannon
## --printshellcmds --verbose

# snakemake --sdm apptainer --apptainer-args --nv --conda-prefix /n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs --snakefile workflow/Snakefile --configfile config/config.yaml --use-conda --workflow-profile profiles/slurm --profile cannon

snakemake --conda-prefix /n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs --snakefile workflow/Snakefile --configfile config/config.yaml --use-conda --workflow-profile profiles/slurm --profile cannon
