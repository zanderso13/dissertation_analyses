function run_subject_firstlevel_dissertation(PID)
%% var set up
if nargin==0 % defaults just for testing 
    % Define some 
    PID = "21684"; 
    
end

are_you_doing_activation_first_levels = 0;

are_you_doing_mid_beta_series = 1;

are_you_doing_rest = 0;

overwrite = 1;
ses = 2;
run = 2;
ndummies = 2;

contrast = 'anticipation'; % consumption

% Define some paths
basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging';

preproc_dir = fullfile(basedir,'smoothed_data/t1/');

numPID = num2str(PID);
PID = strcat('sub-',numPID);


fprintf(['Preparing 1st level model for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


if are_you_doing_mid_beta_series == 1
    timing_dir = fullfile(basedir,'final_betaseries_timing/');
    fl_save_dir = fullfile(basedir,'beta_series/');

    % FL directory for saving 1st level results: beta images, SPM.mat, etc.
    in{1} = {fullfile(fl_save_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)))};

    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-mid_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end

    % onset files
    in{3} = filenames(fullfile(timing_dir, strcat(numPID,'_run-',num2str(run),'.mat')));
    
    if isempty(in{3})
        warning('No modeling found (behav data might be missing)')
        return
    end
    
    in{4} = filenames(fullfile(basedir, 'final_confound/', strcat(numPID,'*mid_run-',num2str(run),'*.mat')));
    
    % check for SPM.mat and overwrite if needed
    skip = 0;
    if exist(fullfile(in{1}{1},'SPM.mat'),'file')
        if overwrite
            fprintf('\n\nWARNING: EXISTING SPM.MATAND BETA FILES WILL BE OVERWRITTEN\n%s\n\n',fullfile(in{1}{1},'SPM.mat'));
            rmdir(in{1}{1},'s');
        else
            fprintf('\n\nFirst levels already exist, wont ovewrite: %s\n\n',fullfile(in{1}{1},'SPM.mat'));
            skip=1;
        end
    end
    
    if ~skip
        % make dir for beta and contrast files
        if ~isdir(in{1}{1}), mkdir(in{1}{1}); end    

        % run spm FL estimation
        cwd = pwd;
        job = strcat('MID_beta_series_template.m');
        %%
        spm('defaults', 'FMRI')
        spm_jobman('serial',job,'',in{:});

        cd(cwd);
    end
    
end

if are_you_doing_rest == 1
    
    ndummies = 10;

    % FL directory for saving 1st level results: beta images, SPM.mat, contrasts, etc.
    in{1} = {fullfile(basedir, 'resting_state_fl', PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)))};

    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end
    
    %% nuisance covs

    % fmriprep output
    confound_fname = filenames(fullfile(basedir, 'final_confound/', strcat(numPID,'*rest_run-1.mat')));
   
    in{3} = {confound_fname{1}};

    % checks
    if any(cellfun( @(x) isempty(x{1}), in))
        in
        error('Some input to the model is missing')
    end
    

    % check for SPM.mat and overwrite if needed
    
    if exist(fullfile(in{1}{1},'SPM.mat'),'file')
        if overwrite
            fprintf('\n\nWARNING: EXISTING SPM.MATAND BETA FILES WILL BE OVERWRITTEN\n%s\n\n',fullfile(in{1}{1},'SPM.mat'));
            rmdir(in{1}{1},'s');
        else
            fprintf('\n\nFirst levels already exist, wont ovewrite: %s\n\n',fullfile(in{1}{1},'SPM.mat'));
            skip=1;
        end
    end


    % run spm FL estimation
    cwd = pwd;
    job = strcat('SPM_rest_template.m');
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

end

%% Model for MID task. First pass at first levels --> activation
if are_you_doing_activation_first_levels == 1
    % FL directory for saving 1st level results: beta images, SPM.mat, contrasts, etc.
    in{1} = {fullfile(basedir, '/activation/', PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)))};
    
    contrast = 'anticipation'; 
    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',num2str(numPID),'.*task-mid_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end
    
    % onset files
    in{3} = filenames(fullfile(basedir,'/mid_spm_timing_baseline/', strcat('run-',num2str(run)), contrast, strcat('*',numPID,'*.mat')));

    if isempty(in{3})
        warning('No modeling found (behav data might be missing)')
        return
    end
    %% nuisance covs

    % fmriprep output
    confound_fname = filenames(fullfile(basedir, 'final_confound/', strcat(numPID,'*mid_run-',num2str(run),'*.mat')));
    
    in{4} = {confound_fname{1}};

    % checks
    if any(cellfun( @(x) isempty(x{1}), in))
        in
        error('Some input to the model is missing')
    end
    
    % check for SPM.mat and overwrite if needed
    
    if exist(fullfile(in{1}{1},'SPM.mat'),'file')
        if overwrite
            fprintf('\n\nWARNING: EXISTING SPM.MATAND BETA FILES WILL BE OVERWRITTEN\n%s\n\n',fullfile(in{1}{1},'SPM.mat'));
            rmdir(in{1}{1},'s');
        else
            fprintf('\n\nFirst levels already exist, wont ovewrite: %s\n\n',fullfile(in{1}{1},'SPM.mat'));
            skip=1;
        end
    end


    % run spm FL estimation
    cwd = pwd;
    job = strcat('SPM_',contrast,'_template.m');
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});
    
end
    
clear in
