doing_beta_series = 1;

if doing_beta_series == 1
    atl = fmri_data('~/repo/Schizconnect/AAL3/AAL3v1.nii');
    vs = fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds/VS_8mmsphere_Oldham_Rew.nii');
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/beta_series/';
    cd(basedir)
    spm_fnames = filenames(fullfile(strcat('sub-*/ses-2/run-1/SPM.mat')));
    for f = 1:length(spm_fnames)
        pid = spm_fnames{f}(5:9);
        fprintf(strcat(pid,'\n'))
        load(spm_fnames{f})
        beta_file_nums = find(contains(SPM.xX.name(:),'antgain'));
        temp_num_strings = {};
        beta_fnames = {};
        for nums = 1:length(beta_file_nums)
            temp_num_strings{nums} = pad(num2str(beta_file_nums(nums)),4,'left','0');
            beta_fnames{nums} = strcat(basedir,'/sub-',pid,'/ses-2/run-1/','beta_',temp_num_strings{nums},'.nii');
        end
        dat = fmri_data(beta_fnames);
        aaldat = extract_roi_averages(dat, atl);
        vsdat = extract_roi_averages(dat, vs);
        
        % probably best to VS x voxel, VS voxel x voxel, VS voxel x AAL
        % regions

        for voxel = 1:length(dat.dat)
            corr_vs_to_wholebrain_voxel(voxel,:) = corr(dat.dat(voxel,:)',vsdat.dat);
            
        end

        for region = 1:length(aaldat)
            corr_vs_to_wholebrain_aal(region,1) = corr(aaldat(region).dat,vsdat.dat);
        end

        for region = 1:length(aaldat)
            for voxel = 1:length(vsdat.all_data)
                corr_vs_voxel_to_region(region,voxel) = corr(aaldat(region).dat,vsdat.all_data(:,voxel));
            end
        end
        
        curr_filename = strcat(pid,'_beta_correlations_and_extracted_mid_data.mat');

        save(fullfile(basedir,curr_filename),"corr_vs_voxel_to_region", "corr_vs_to_wholebrain_aal", "corr_vs_to_wholebrain_voxel", "dat")

    end

end