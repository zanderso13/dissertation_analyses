function prep_connectivity_matrices(pid,run)

doing_beta_series = 0;
doing_rest = 1;

if doing_beta_series == 1
    atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    % atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');
    vs = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/VS_8mmsphere_Oldham_Rew.nii');
    ofc = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/OFC_8mmsphere_Oldham.nii');
    acc = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/bilat_AAL3_ACC.nii');
    amyg = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/HO_Amygdala_50prob.nii');
    
    savedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices';
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/beta_series';
    cd(basedir)
%     spm_fnames = filenames(fullfile(strcat('sub-*/ses-2/gain_contrasts/run-2/SPM.mat')));
%     keyboard
%     for f = 1:length(spm_fnames)

%         pid = spm_fnames{f}(5:9);
        fprintf(['Preparing rest corr for ' num2str(pid)]);

        curr_file = filenames(fullfile(strcat('*',num2str(pid),'/'),'*/gain_contrasts/',run,'SPM.mat'));
        load(curr_file{1});
        beta_file_nums = find(contains(SPM.xX.name(:),'antgain'));
        temp_num_strings = {};
        beta_fnames = {};
        for nums = 1:length(beta_file_nums)
            temp_num_strings{nums} = pad(num2str(beta_file_nums(nums)),4,'left','0');
            beta_fnames{nums} = strcat(basedir,'/sub-',num2str(pid),'/ses-2/gain_contrasts/',run,'/beta_',temp_num_strings{nums},'.nii');
        end
        dat = fmri_data(beta_fnames);
        seitz = extract_roi_averages(dat, atl);
        vsdat = extract_roi_averages(dat, vs);
        accdat = extract_roi_averages(dat, acc);
        ofcdat = extract_roi_averages(dat, ofc);
        amygdat = extract_roi_averages(dat, amyg);
        
        % first for loop gives seed to whole brain as a sanity check for
        % different networks related to the 4 seeds I have

        for voxel = 1:length(dat.dat)
            corr_vs_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',vsdat.dat);  
            corr_ofc_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',ofcdat.dat); 
            corr_acc_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',accdat.dat); 
            corr_amyg_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',amygdat.dat); 
        end
        
        % next for loop creates a 300x300 corr matrix
        for r = 1:length(seitz)
            seitz_mat(:,r) = seitz(r).dat;
        end
       

        % final loops need to extract seed to 300 regions corr mats
        for r = 1:size(seitz_mat,2)
            for voxel = 1:length(vsdat.all_data)
                corr_vs_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),vsdat.all_data(:,voxel));
            end
        end

        for r = 1:size(seitz_mat,2)
            for voxel = 1:length(ofcdat.all_data)
                corr_ofc_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),ofcdat.all_data(:,voxel));
            end
        end

        for r = 1:size(seitz_mat,2)
            for voxel = 1:length(accdat.all_data)
                corr_acc_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),accdat.all_data(:,voxel));
            end
        end

        for r = 1:size(seitz_mat,2)
            for voxel = 1:length(amygdat.all_data)
                corr_amyg_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),amygdat.all_data(:,voxel));
            end
        end
        
        curr_filename = strcat(num2str(pid),'_',run,'.mat');

        save(fullfile(savedir,curr_filename),"corr_vs_wholebrain", "corr_ofc_wholebrain",...
            "corr_amyg_wholebrain", "corr_acc_wholebrain","seitz_mat",...
            "corr_vs_voxel_to_region","corr_ofc_voxel_to_region","corr_acc_voxel_to_region",...
            "corr_amyg_voxel_to_region")

    %end

end

if doing_rest == 1

    atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    % atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');
    vs = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/VS_8mmsphere_Oldham_Rew.nii');
    ofc = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/OFC_8mmsphere_Oldham.nii');
    acc = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/bilat_AAL3_ACC.nii');
    amyg = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/HO_Amygdala_50prob.nii');
    
    savedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices';
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/AIB_RestOnly_Newtworks';
    cd(basedir)
    
    rest_fname = filenames(fullfile(strcat('sub-',num2str(pid),'*/ses-1/ssub*.nii')));

    dat = fmri_data(rest_fname{1});
    seitz = extract_roi_averages(dat, atl);
    vsdat = extract_roi_averages(dat, vs);
    accdat = extract_roi_averages(dat, acc);
    ofcdat = extract_roi_averages(dat, ofc);
    amygdat = extract_roi_averages(dat, amyg);
    

    % first for loop gives seed to whole brain as a sanity check for
    % different networks related to the 4 seeds I have

    for voxel = 1:length(dat.dat)
        corr_vs_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',vsdat.dat);  
        corr_ofc_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',ofcdat.dat); 
        corr_acc_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',accdat.dat); 
        corr_amyg_wholebrain(voxel,:) = corr(dat.dat(voxel,:)',amygdat.dat); 
    end
    
    % next for loop creates a 300x300 corr matrix
    for r = 1:length(seitz)
        seitz_mat(:,r) = seitz(r).dat;
    end
   
    
    % final loops need to extract seed to 300 regions corr mats
    for r = 1:size(seitz_mat,2)
        for voxel = 1:size(vsdat.all_data,2)
            corr_vs_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),vsdat.all_data(:,voxel));
        end
    end

    for r = 1:size(seitz_mat,2)
        for voxel = 1:size(ofcdat.all_data,2)
            corr_ofc_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),ofcdat.all_data(:,voxel));
        end
    end

    for r = 1:size(seitz_mat,2)
        for voxel = 1:size(accdat.all_data,2)
            corr_acc_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),accdat.all_data(:,voxel));
        end
    end

    for r = 1:size(seitz_mat,2)
        for voxel = 1:size(amygdat.all_data,2)
            corr_amyg_voxel_to_region(r,voxel) = corr(seitz_mat(:,r),amygdat.all_data(:,voxel));
        end
    end
    
    curr_filename = strcat(num2str(pid),'_rest.mat');
    how_many_volumes = size(dat.dat,2);
    save(fullfile(savedir,curr_filename),"corr_vs_wholebrain", "corr_ofc_wholebrain",...
        "corr_amyg_wholebrain", "corr_acc_wholebrain","seitz_mat",...
        "corr_vs_voxel_to_region","corr_ofc_voxel_to_region","corr_acc_voxel_to_region",...
        "corr_amyg_voxel_to_region", "how_many_volumes")
    

end
