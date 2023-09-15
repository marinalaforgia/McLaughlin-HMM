#!/bin/bash -l
#
#SBATCH -n 12 #number of cores
#SBATCH -D simulation/ #location on cluster of files you are using
#SBATCH -J McL-hmm-sim #name of job
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=mlaforgia@ucdavis.edu # send-to address
#SBATCH --mem-per-cpu=2G #memory per node in Gb
#SBATCH -t 48:00:00 #time in hours:min:sec

module load R

R CMD BATCH McLaughlin_simulation.R



