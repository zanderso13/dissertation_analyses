fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices/combined_runs_with_hyp/*mat'));

% only one of these can be turned 'on' at a time. Composites should be
% used when you're 1) trying to optimize speed or 2) want to get voxelwise
% translation. Basically, this feeds hyperalignment matrices that yield a
% 162x163 rotation matrix, the less computationally expensive route. On the
% other hand, if you need voxelwise rotation information, then grab
% rotation. It will take longer and it will be important to change the perm
% loop to be something like 1:10. 1:1000 will crash your computer
composites = 0;
rotation = 1;
% to optimize signal for plotting, I'm going to feed hyperalignment
% controls first and psychosis second. I'll randomize the order of subjects
% across 10 permutations, but am trying to stack the deck to build a
% hyperaligned signature of psychosis in this sample
pure_random = 1;
within_diagnosis = 0;
% when rotation = 1, this value should be small (like 10). The rotation
% method mentioned above is super expensive computationally. If you're
% going for composites then feel free to set this up at 1000. That
% shouldn't crash your computer and you'll want the extra iterations. This
% is especially true if you're doing a pure random test, which will shuffle
% subjects in all kinds of different orders.
permutation_length = 100;

load(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_motion_exclusions.mat'));

% apply exclusions based on a >0.2mm FD 
for exclude = 1:length(pid_exclude_list)
    fnames(contains(fnames(:),pid_exclude_list{exclude,1})) = [];
end

for f = 1:length(fnames)
    load(fnames{f});
    unaligned_mats{f} = hypmat;
end
unaligned_mats2 = unaligned_mats;

fprintf(strcat('Im planning to do this many permutations:',num2str(permutation_length),'\n'))

for perm = 1:permutation_length
    if pure_random == 1
        idx = randperm(length(unaligned_mats2));
    
        unaligned_data=unaligned_mats2(idx);

        [aligned_data,transforms] = hyperalign(unaligned_data{:});
    
        [~,sidx]=sort(idx);
    
        transforms_resorted=transforms(sidx);    
        aligned_data_resorted=aligned_data(sidx);
    elseif within_diagnosis == 1
        clin_data = unaligned_mats2(dat.Diagnosis==1);
        hc_data = unaligned_mats2(dat.Diagnosis==0);
        
        idx_clin = randperm(length(clin_data));
        idx_hc = randperm(length(hc_data));

        data = [hc_data(idx_hc),clin_data(idx_clin)];
        [aligned_data,transforms] = hyperalign(data{:});

        [~,sidx_clin]=sort(idx_clin);
        [~,sidx_hc]=sort(idx_hc);

        transforms_resorted_clin = transforms(length(hc_data)+1:length(transforms)); 
        transforms_resorted_clin = transforms_resorted_clin(sidx_clin);

        transforms_resorted_hc = transforms(1:length(hc_data));
        transforms_resorted_hc = transforms_resorted_hc(sidx_hc);

        transforms_resorted = transforms_resorted_clin;
    end
    fprintf(strcat('Im on permutation: ',num2str(perm),'\n'))
%     curr_filename = fullfile(strcat('/Users/zacharyanderson/Documents/ADAPTlab/schizconnect/final_data/all_permutations_aligned_data/aligned_data_perm',num2str(perm),'.mat'));
%     save(curr_filename, 'aligned_data_resorted');
    % this bit calculates and stores composite values
    for sub = 1:length(transforms_resorted)
        % magnitude of translation (distance metric)
        translation_temp = transforms_resorted{sub}.c.^2;
        translation(sub,perm)= sqrt(sum(translation_temp(:)));
        translation_all{sub}(perm,:) = sum(abs(transforms_resorted{sub}.c),2);
%         
%         % magnitude of rotation (distance metric)
        rotation_temp = transforms_resorted{sub}.T.^2;
        rotation(sub,perm) = sqrt(mean(rotation_temp(:)));
        rotation_all{sub}(perm,:) = sum(abs(transforms_resorted{sub}.T),2);
%         % magnitude of scaling (single value)
        scale(sub,perm) = transforms_resorted{sub}.b;
%         % determinant of rotation matrix (were data rotated with
%         % reflection (det(T) = -1, or not?
        reflection(sub,perm) = det(transforms_resorted{sub}.T);
%         frontal_vox_temp = transforms_resorted{sub}.c * transforms_resorted{sub}.T';
%         frontal_translation(sub,:) = frontal_vox_temp(1,:)';
%         frontal_rot(sub,:) = mean(transforms_resorted{sub}.T,2)';
        
    end
    %curr_filename = fullfile('/Users/zacharyanderson/Documents/ADAPTlab/schizconnect/final_data/',strcat('rotation_single_perm_',num2str(perm),'.mat'));
    %save(curr_filename,'frontal_translation','frontal_rot','-v7.3')
end

save hyperalignment_composite_perm.mat translation rotation scale reflection translation_all rotation_all -v7.3