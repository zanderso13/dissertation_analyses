#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH --mem=64G

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('~/repo')); gica_cmd --data ~/repo/subject_list_for_ica.txt --o '/projects/p30954/brainmapd_ica'; quit"

