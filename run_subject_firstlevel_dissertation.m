function run_subject_firstlevel_BrainMAPD_PPI(PID)
%% var set up
if nargin==0 % defaults just for testing 
    % Define some 
    PID = "21684"; 
    
end

are_you_doing_activation_first_levels = 0;
% if you're doing resting state analysis, next line should be set to 0
are_you_doing_ppi_first_levels = 0;
% Rod suggested I do traditional resting state analysis on task to show
% that reward networks are active. This will also work better for the way I
% do hyperalignment
are_you_doing_denoise_mid = 1;

% the next two lines are a bit redundant, but it's how I've tried to get
% around the different naming conventions that I'm used to seeing. There
% are places throughout this script that reference these strings as part of
% file names that SPM is either reading in or outputting. 

task = 'mid'; % 'rest', 'mid'
contrast = 'anticipation'; % anticipation, rest, consumption

% the next line only applies if you're doing ppi
seed_region = 'Oldham_Rew'; % anticipation: Amygdala, OFC, Oldham_Rew (VS), Oldham_Loss (VS); consumption: Amygdala, OFC, Oldham_Con (VS)
overwrite = 1;
ses = 2;
run = 1;
ndummies = 2; % 10 for rest, 2 for mid

% Define some paths
basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/zach_and_nina_first_levels/';

% directories
% first is where your activation related stats files will be output to. For
% rest, change it to rest! For mid change it to activation. For mid
% denoising/traditional resting state analysis on task, mid_denoise
fl_dir = fullfile(basedir,'/mid_denoise');
% next is where the preprocessed data is
preproc_dir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/smoothed_functional_data';
% where framewise displacement files will be saved
save_dir = fullfile(basedir,'/FD');
% directory where I'm storing timing files for the MID
timing_dir = fullfile(strcat('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_spm_timing_baseline/run-',num2str(run)),contrast);
% this is where the ppi specific models will be output
ppi_fl_dir = fullfile(basedir,'/ppi');
% this is where masks for the current study are held
seed_dir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds';
% this is where confound files are. these are distilled separately and then
% implemented here
confound_dir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/first_level_confounds_dissertation';

numPID = num2str(PID);
PID = strcat('sub-',numPID);


fprintf(['Preparing 1st level model for ' task ' task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);



%% Model for MID task. First pass at first levels --> activation
if are_you_doing_activation_first_levels == 1
    % FL directory for saving 1st level results: beta images, SPM.mat, contrasts, etc.
    in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)))};

    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-', task,'_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end
    
    if strcmp(task,'mid') == 1
        % onset files
        in{3} = filenames(fullfile(timing_dir, strcat(numPID,'*')));
    
        if isempty(in{3})
            warning('No modeling found (behav data might be missing)')
            return
        end
        %% nuisance covs
    
        % fmriprep output
        confound_fname = filenames(fullfile(confound_dir, strcat(PID,'*',task,'_ses-2_run-',num2str(run),'*confounds*.mat')));
	    
        % choose which matrix to use
        load(confound_fname{1});
    
       
        in{4} = {confound_fname{1}};
    
        % checks
        if any(cellfun( @(x) isempty(x{1}), in))
            in
            error('Some input to the model is missing')
        end
    else
        %% nuisance covs
    
        % fmriprep output
        confound_fname = filenames(fullfile(confound_dir, strcat(PID,'*',task,'_ses-2_run-',num2str(run),'*confounds*.mat')));
    	
        % choose which matrix to use
        load(confound_fname{1});
    
       
        in{3} = {confound_fname{1}};
    
        % checks
        if any(cellfun( @(x) isempty(x{1}), in))
            in
            error('Some input to the model is missing')
        end
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
    
    
if are_you_doing_ppi_first_levels == 1
    % FL directory for saving 1st level results: beta images, SPM.mat, etc.
    in{1} = {fullfile(ppi_fl_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)))};

    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-mid_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end

    % onset files
    in{3} = filenames(fullfile(timing_dir, strcat(numPID,'*')));

    if isempty(in{3})
        warning('No modeling found (behav data might be missing)')
        return
    end
    
    regressor_temp = load(fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)),'SPM.mat')); 
    % choose which matrix to use
    R = regressor_temp.SPM.xX.X(:,[1:4,6:size(regressor_temp.SPM.xX.X,2)-1]);
    % put in SPM format: matrix called 'R', and 'names'
    names = regressor_temp.SPM.xX.name(:,[1:4,6:size(regressor_temp.SPM.xX.X,2)-1]);

    % use canlab core tools to pull out ROI time course data/convolve
    dat = fmri_data(filenames(fullfile(preproc_dir, strcat('ssub-',numPID,'*mid*run-',num2str(run),'*'))));
    roi = fmri_data(filenames(fullfile(seed_dir,strcat('*',seed_region,'*nii'))));
    roi_dat = extract_roi_averages(dat,roi);
    
    % multiply time course by each original task regressor. These are the
    % ppi regressors that we will analyze in the final model

    R(:,1) = R(:,1) .* zscore(roi_dat.dat(ndummies+1:size(roi_dat.dat,1))); names{1} = 'lose_ppi';
    R(:,2) = R(:,2) .* zscore(roi_dat.dat(ndummies+1:size(roi_dat.dat,1))); names{2} = 'lose_0_ppi';
    R(:,3) = R(:,3) .* zscore(roi_dat.dat(ndummies+1:size(roi_dat.dat,1))); names{3} = 'gain_ppi';
    R(:,4) = R(:,4) .* zscore(roi_dat.dat(ndummies+1:size(roi_dat.dat,1))); names{4} = 'gain_0_ppi';

    % include the average seed time course in the model 
    R = [R,zscore(roi_dat.dat(ndummies+1:size(roi_dat.dat,1)))]; names{length(names)+1} = 'timecourse';   
    
    confound_fname_ppi = fullfile(confound_dir, strcat(PID,'_mid_ses-2_run-',num2str(run),'_ppi_confounds.mat'));
    
    save(confound_fname_ppi, 'R','names')
    
    in{4} = {confound_fname_ppi};

    % checks
    if any(cellfun( @(x) isempty(x{1}), in))
        in
        error('Some input to the model is missing')
    end
    
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

        save(fullfile(ppi_fl_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)),'mid_confounds.mat'),'R','names');


        % run spm FL estimation
        cwd = pwd;
        job = strcat('MID_SPM_PPI_',contrast,'_template.m');
        %%
        spm('defaults', 'FMRI')
        spm_jobman('serial',job,'',in{:});

        cd(cwd);
    end
end

if are_you_doing_denoise_mid == 1

    % FL directory for saving 1st level results: beta images, SPM.mat, contrasts, etc.
    in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), contrast, strcat('run-', num2str(run)))};

    % preproc images
    in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-', task,'_run-',num2str(run),'_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii'), ndummies+1:9999));

    if isempty(in{2}{1})
        warning('No preprocd functional found')
        return
    end
    
    %% nuisance covs

    % fmriprep output
    confound_fname = filenames(fullfile(confound_dir, strcat(PID,'*',task,'_ses-2_run-',num2str(run),'*confounds*.mat')));
	
    % choose which matrix to use
    load(confound_fname{1});

   
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
    job = strcat('SPM_mid_denoise_template.m');
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

end


