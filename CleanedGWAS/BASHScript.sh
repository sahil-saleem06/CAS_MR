#!/bin/bash
#BSUB -J Clean_PreLift_GWAS
#BSUB -n 1
#BSUB -R rusage[mem=64000]
#BSUB -q voltron_normal
#BSUB -W 48:00
#BSUB -o Clean_PreLift_GWAS.out
#BSUB -e Clean_PreLift_GWAS.err

# Set R library directory
export R_LIBS_USER=$HOME/R/rocker-rstudio/bioconductor-tidyverse_3.17

# Ensure singularity sees the correct folders
export SINGULARITY_BIND="/lsf:/lsf, /project/:/project/, /appl/:/appl/, /scratch/:/scratch, /static:/static"

# Stop script on first error
set -e

# Load singularity module
module load singularity

# Execute the R script
singularity exec \
  /project/voltron/rstudio/containers/bioconductor-tidyverse_3.17.sif \
  Rscript /project/damrauer_scratch/Users/saleemsa/CleanedGWAS/RScriptCarotid.R
