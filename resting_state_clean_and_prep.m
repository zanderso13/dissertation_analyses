% this script is intended to be run directly after first levels to prepare
% resting state data for group analysis

basedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS';
atlasdir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/AAL3';

cd(basedir)

fnames = filenames(fullfile('sub*/ses-2/rest/run-1/SPM.mat'));
atl = fmri_data(filenames(fullfile(atlasdir,'AAL3v1.nii')));

% apply AAL3
for sub = 1:length(fnames)
    pid = fnames{sub}(5:9);
    dat = fmri_data(filenames(fullfile(strcat('*',pid,'*/ses-2/rest/run-1/Res*'))));
    r{sub} = extract_roi_averages(dat,atl);
end

% create basic connectivity matrix
for sub = 1:length(r)
    for region = 1:length(r{sub})
        conn{sub}(:,region) = r{sub}(region).dat(:);
    end

    corr_matrices{sub} = corr(conn{sub});
end





