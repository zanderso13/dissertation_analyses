#!/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 08:00:00
#SBATCH --mem=30G

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('/home/zaz3744/repo')); dissertation_main_analysis_gendis_rest; quit"

