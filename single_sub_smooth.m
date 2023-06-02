
function smooth_single_sub(PID,ses,run,overwrite)
%% var set up
if nargin==0 % defaults just for testing
    PID = 10001;  
    overwrite = 1;
    ses = 2;
    run = 1;
end

preproc_dir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/fmriprep';

if nargin==1
    overwrite = 1;
end  

PID = strcat('sub-',num2str(PID));
ndummies=0;
% FL directory for saving 1st level results: beta images, SPM.mat, etc.
% in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};

rundir = fullfile(preproc_dir, PID, strcat('ses-', num2str(ses)), 'func');

in{1} = cellstr(spm_select('ExtFPList', rundir, strcat('.*task-REST_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(in{1}{1})
    warning('No preprocd functional found')
    return
end

jobfile = {'/home/zaz3744/repo/acnlab_repo/preproc/smooth_template.m'};
jobs = 'smooth_template.m';


spm('defaults', 'FMRI');
spm_jobman('run', jobs, in{:});


end
