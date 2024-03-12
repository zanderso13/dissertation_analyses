basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging';
scriptdir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/scriptdir';
corrdir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices';
repodir = '/home/zaz3744/repo';

% What run of your task are you looking at?
run = 'run-1'; % run-2, rest will only have run-1

% What task are you running?
task = 'AIB_RestOnly_Newtworks'; % eventually beta_series 

cd(fullfile(basedir,task))

fnames = filenames(fullfile('sub*/ses-1/ssub*rest_final.nii'));

conn_list = filenames(fullfile(corrdir,strcat('*.mat')));

counter = 1;
for sub = 1:length(fnames)
    
    curr_sub = fnames{sub}(5:9);
    if isempty(find(contains(conn_list,curr_sub)))
        new_list(counter) = str2num(curr_sub);
        counter = counter + 1;
    else
        continue
    end
end

% loop through the new list of subjects and submit each prep connectivity
% job separately

cd(scriptdir)
keyboard
for sub = 1:length(new_list)
    ids = new_list(sub);

    % this bit will run the loop one subject at a time for testing

%     prep_connectivity_matrices_2(num2str(ids),run)

%     % this bit will write out text files then submit them to the cluster
        s = ['#!/bin/bash\n\n'...
     '#SBATCH -A p30954\n'...
     '#SBATCH -p short\n'...
     '#SBATCH -t 01:00:00\n'...  
     '#SBATCH --mem=30G\n\n'...
     'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); prep_connectivity_matrices_2(' num2str(ids) ',''' run '''); quit"\n\n']; 
   
     scriptfile = fullfile(scriptdir, 'first_level_script.sh');
     fout = fopen(scriptfile, 'w');
     fprintf(fout, s);
     
     !chmod 777 first_level_script.sh
     !sbatch first_level_script.sh
%      

end
