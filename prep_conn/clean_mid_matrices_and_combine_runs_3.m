% I need to combine runs and apply exclusions. I will output single files
% that contain both runs of the MID matrices. 
basedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/mid_corr_matrices';
contrast = 'antgain'; % hitcongain
atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');

hyper = 0;

cd(basedir)

fnames1 = filenames(fullfile('*_run-1*'));
fnames2 = filenames(fullfile('*_run-2*'));
load(fnames2{1});
dat = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/sub-10001/ses-2/gain_contrasts/run-1/beta_0001.nii');

for f = 1:length(fnames2)
    pid = fnames2{f}(1:5);
    fprintf(strcat('Working on: ',pid, '\n\n\n'))
    if sum(contains(fnames1(:),pid)) ==1
        % first, let's just look at all the brains I have both runs for
        m1.(strcat('sub',pid)) = load(fnames1{contains(fnames1(:),pid)});
        m2.(strcat('sub',pid)) = load(fnames2{f});       

        % this is where you can grab different seeds. here are the examples
        % you can use: corr_acc_wholebrain corr_ofc_wholebrain
        % corr_vs_wholebrain corr_amyg_wholebrain
        temp_brain(:,f) = m1.(strcat('sub',pid)).corr_vs_wholebrain; 
        temp_brain2(:,f) = m2.(strcat('sub',pid)).corr_vs_wholebrain;
        isnan(temp_brain(:,f)) = 0;
        isnan(temp_brain2(:,f)) = 0;
    end
end

dat.dat = temp_brain; dat.X = ones(size(temp_brain,2),1);

% average brain, then threshold at 0.2

dat.dat = mean(dat.dat,2); dat.dat(dat.dat<0.2)=0;

