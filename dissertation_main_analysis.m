analyze_with_behavioral_and_self_report = 1;
hyper = 1;
%% below this is old code I need to repurpose for the new sets of analyses I'm running. 
% beta series correlations will be used for task. resting state func conn
% will be what it is. hyperalignment section will be appended at the very
% end. 
basedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices/combined_runs_with_hyp';
cd(basedir)
if analyze_with_behavioral_and_self_report == 1
    fnames = filenames('*.mat');

    load(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_motion_exclusions.mat'));

    % apply exclusions based on a >0.2mm FD 
    for exclude = 1:length(pid_exclude_list)
        fnames(contains(fnames(:),pid_exclude_list{exclude,1})) = [];
    end


    %% prep behavioral data and match to pids
    % behavioral data load in
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/demographics.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/trilevel.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/immune.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/meds.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/chronic_lsi_t2.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/chronic_lsi_t1.mat');
    load('/Users/zacharyanderson/Documents/dissertation/outcomes/income.mat');
    
    
    for sub = 1:length(fnames)
        id{sub,1} = str2num(fnames{sub}(1:5));
        % immune
        if ~isempty(find(immune.PID==id{sub,1}))
            i(sub,1) = immune.T1BDicsavg(find(immune.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing immune data \n'))
            i(sub,1) = NaN;
        end
        % t1 lsi
        if ~isempty(find(BrainMAPDT1LSIChronic202047.PID==id{sub,1}))
            traum1(sub,1) = sum(table2array(BrainMAPDT1LSIChronic202047(find(BrainMAPDT1LSIChronic202047.PID==id{sub,1}),6:19)));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing immune data \n'))
            traum1(sub,1) = NaN;
        end
        % t2 lsi
        if ~isempty(find(BrainMAPDT2LSIChronic2020221.PID==id{sub,1}))
            traum2(sub,1) = sum(table2array(BrainMAPDT2LSIChronic2020221(find(BrainMAPDT2LSIChronic2020221.PID==id{sub,1}),6:19)));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing traum1 data \n'))
            traum2(sub,1) = NaN;
        end
        % income
        if ~isempty(find(BrainMAPDT1AdultSES2021129.PID==id{sub,1}))
            inc(sub,1) = BrainMAPDT1AdultSES2021129.T1SESfaminc(find(BrainMAPDT1AdultSES2021129.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing traum2 data \n'))
            inc(sub,1) = NaN;
        end
        % demographics
        if ~isempty(find(dem.PID==id{sub,1}))
            d(sub,:) = dem(find(dem.PID==id{sub,1}),:);
        else
            fprintf(strcat(num2str(id{sub,1}),': missing dem data \n'))
            d(sub,1) = NaN;
        end
        % meds
        if ~isempty(find(meds.PID==id{sub,1}))
            m(sub,1) = meds.T2SCcmidopany(find(meds.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing med data \n'))
            m(sub,1) = 0;
        end
        % symptoms
        if ~isempty(find(trilevel.PID==id{sub,1}))
            s(sub,1) = trilevel.T2gendis(find(trilevel.PID==id{sub,1}));
            s(sub,2) = trilevel.T2anhedon(find(trilevel.PID==id{sub,1}));
            s(sub,3) = trilevel.T2fears(find(trilevel.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing clinical data \n'))
            s(sub,1) = NaN;
            s(sub,2) = NaN;
            s(sub,3) = NaN;
        end
        
    end
    % site, age, sex, racial and ethnic identity, and participant reports
    % of annual family income - covariates from greg's paper.
    site = cell2mat(id); site(site<20000) = 0; site(site>20000) = 1; site = array2table(site);
    s = array2table(s); s.Properties.VariableNames = {'GeneralDistress','AnhedoniaApprehension','Fears'};
    m = array2table(m); m.Properties.VariableNames = {'medication'};
    i = array2table(i); i.Properties.VariableNames = {'inflammation'};
    %d = array2table(d); d.Properties.VariableNames = dem.Properties.VariableNames;
    
    final_outcomes = [s,m,i,d,site,array2table(traum1),array2table(traum2),array2table(inc)];
    
    regressors_variablenames = {'inflammation','sex','site','medication','ethnicity','race','lsi1','lsi2','famincome'};
    regressors = [final_outcomes.inflammation,final_outcomes.sex,final_outcomes.site,final_outcomes.medication,final_outcomes.ethnicity,final_outcomes.race];
    
    %% prep brain data
    % the sections commented out below were run once and then I saved the
    % final data object so I wouldn't have to regenerate it. It's still
    % there so I can go back and recreate if needed.

    load(fnames{1});
    final_brain = temp_brain;
    for f = 1:length(fnames)
        load(fnames{f});
        if sum(isnan(corr_vs_to_wholebrain_voxel)) > 0
            corr_vs_to_wholebrain_voxel(isnan(corr_vs_to_wholebrain_voxel)) = 0;
            final_brain.dat(:,f) = corr_vs_to_wholebrain_voxel;
        else
            final_brain.dat(:,f) = corr_vs_to_wholebrain_voxel;
        end
    end
    %load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_vs_to_connectivity.mat');
    final_brain.X = regressors;
    regressor_ind = ~isnan(final_brain.X(:,1));
    final_brain.dat(:, isnan(final_brain.X(:,1))) = [];
    final_brain.X(isnan(final_brain.X(:,1)),:) = [];
    statobj = regress(final_brain);
    threshobj = threshold(statobj.t,0.005,'unc','k',10);

    %% seed to seed
    atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');
    % atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/AAL3/AAL3v1.nii');
    vs_to_region = extract_roi_averages(final_brain,atl);

    for a = 1:length(vs_to_region)
        temp = fitlm(final_brain.X, vs_to_region(a).dat);
        if temp.Coefficients.pValue(2) < 0.05
            a
        end

    end    


    if hyper == 1
        
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
        
        fnames2 = fnames(regressor_ind);
        for f = 1:length(fnames2)
            load(fnames2{f});
            hypmat(isnan(hypmat)) = 0;
            unaligned_mats{f} = hypmat';
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
            outcome = final_brain.X(:,1);

            for i=1:length(aligned_data_resorted)
                aligned_x(:,i) = aligned_data_resorted{i}(:);
                unaligned_x(:,i) = unaligned_data{i}(:);
            end

              
            [uXL{perm},uYL{perm},uXS{perm},uYS{perm},uBETA{perm},uPCTVAR{perm},uMSE{perm},ustats{perm}] = plsregress(unaligned_x',outcome,20,'CV',10);
            [aXL{perm},aYL{perm},aXS{perm},aYS{perm},aBETA{perm},aPCTVAR{perm},aMSE{perm},astats{perm}] = plsregress(aligned_x',outcome,20,'CV',10);    
            %curr_filename = fullfile('/Users/zacharyanderson/Documents/ADAPTlab/schizconnect/final_data/',strcat('rotation_single_perm_',num2str(perm),'.mat'));
            %save(curr_filename,'frontal_translation','frontal_rot','-v7.3')
        end
        
        save hyperalignment_composite_perm.mat translation rotation scale reflection translation_all rotation_all -v7.3    
        
    end
end


