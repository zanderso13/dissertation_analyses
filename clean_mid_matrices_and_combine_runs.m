% I need to combine runs and apply exclusions. I will output single files
% that contain both runs of the MID matrices. 
basedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices';
contrast = 'antgain'; % hitcongain
atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');

hyper = 0;

cd(basedir)

fnames1 = filenames(fullfile(strcat('*_run-1*',contrast,'*')));
fnames2 = filenames(fullfile(strcat('*_run-2*',contrast,'*')));
load(fnames2{1});


for f = 1:length(fnames2)
    pid = fnames2{f}(1:5);
    fprintf(strcat('Working on: ',pid, '\n\n\n'))
    if sum(contains(fnames1(:),pid)) ==1
        
        m1.(strcat('sub',pid)) = load(fnames1{contains(fnames1(:),pid)});
        m2.(strcat('sub',pid)) = load(fnames2{f});            
        temp_brain = m1.(strcat('sub',pid)).dat;
        temp_brain.dat = [temp_brain.dat,m2.(strcat('sub',pid)).dat.dat];
        % plot(temp_brain)
        
        % outliers = ans;
        %temp_brain.dat = temp_brain.dat(:,~outliers);
        
        vs = fmri_data('/Users/zacharyanderson/Documents/ACNlab/masks/VS_8mmsphere_Oldham_Rew.nii'); % 

        roi = extract_roi_averages(temp_brain,vs);

        other_regions = extract_roi_averages(temp_brain,atl);
        
        for voxel = 1:length(dat.dat)
            corr_vs_to_wholebrain_voxel(voxel,:) = corr(temp_brain.dat(voxel,:)',roi.dat);        
        end
        
        for voxel = 1:length(roi.all_data)
            for reg = 1:length(other_regions)
                hypmat(voxel,reg) = corr(roi.all_data(:,voxel),other_regions(reg).dat); 
            end
        end

        final_save_name = fullfile(strcat('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices/combined_runs_with_hyp/',pid,'_',contrast,'.mat'));

        save(final_save_name,"corr_vs_to_wholebrain_voxel","temp_brain","hypmat");

    end
end


