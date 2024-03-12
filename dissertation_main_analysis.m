% I think this one should always be equal to 1
prep_behavioral_data = 1;

% what task?
mid = 1; rest = 0;

region_name_for_wholebrain_analysis = 'vs'; % vs amyg acc ofc

whole_brain_networks = 0; overwrite_nii = 0; 

linear_seed_to_seed = 1; seeddir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds';
mediation_seed_to_seed = 0; % linear seed to see must be on
moderated_mediation_seed_to_seed = 0;

whole_brain_mediation_analysis = 0; 
whole_brain_moderated_mediation = 0;

hyper= 0; overwrite_hyper = 0; % hardcoded which seed you're looking at

hyper_analyze = 0; remake_multivar_input = 0; do_pls_regress = 0; 

compare_aligned_pls_loadings = 0;

look_at_multivariate_rotation_estimates = 0;

mediate_with_pls_components = 0;

network_based_analyses = 0;

if mid == 1
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats';
end

if rest == 1
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices';
end    

cd(basedir)

fnames = filenames('*.mat');
if mid==1   
    load(fullfile('/home/zaz3744/repo/dissertation_analyses/final_motion_exclusion.mat'));
end

if rest==1
    load(fullfile('/home/zaz3744/repo/dissertation_analyses/final_motion_exclusions_rest.mat'))
end

% apply exclusions based on a >0.2mm FD MID or rest
for exclude = 1:length(pid_exclude_list)
    fnames(contains(fnames(:),pid_exclude_list{exclude,1})) = [];
end

% TASK SPECIFIC EXCLUSIONS
% exclude 1 subject based on NaNs in their connectivity matrix. I think
% this suggests the person didn't have as much variability as they would
% need to generate corrs. This was only true for one region, though I worry
% this reflects a larger issue, potentially related to too much filtering
% of their signal
if mid == 1
    fnames(185) = [];
end

if prep_behavioral_data == 1
    
    %% prep behavioral data and match to pids
    % behavioral data load in
    clinical_base = '/projects/b1108/studies/brainmapd/data/processed/clinical/outcomes';
    
    load(fullfile(clinical_base,'cti.mat'));
    load(fullfile(clinical_base,'trilevel_longitudinal.mat'));
    load(fullfile(clinical_base,'demographics.mat'));
    load(fullfile(clinical_base,'immune_data.mat'));
    load(fullfile(clinical_base,'meds.mat'));
    load(fullfile(clinical_base,'chronic_lsi_t2.mat'));
    load(fullfile(clinical_base,'chronic_lsi_t1.mat'));
    load(fullfile(clinical_base,'income.mat'));
    
    
    for sub = 1:length(fnames)
        id{sub,1} = str2num(fnames{sub}(1:5));
        % immune
        if ~isempty(find(immune.PID==id{sub,1}))
            %i(sub,1) = immune.T1BDicsavg(find(immune.PID==id{sub,1}));
            i(sub,1) = immune.ZT1BDcrpLN(find(immune.PID==id{sub,1})); %
            %p=0.047 for crp-->vs-ofc connectivity in mid
            %i(sub,1) = immune.ZT1BDil8LN(find(immune.PID==id{sub,1})); % p=0.016 in mid vs-acc
            %i(sub,1) = immune.ZT1BDil10LN(find(immune.PID==id{sub,1})); % p=
            
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
        % cti
        if ~isempty(find(cti.subID==id{sub,1}))
            cti_final(sub,:) = cti(find(cti.subID==id{sub,1}),:);
        else
            fprintf(strcat(num2str(id{sub,1}),': missing traum1 data \n'))
            cti_final(sub,:) = array2table(ones(1,size(cti_final,2)));
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
            m(sub,1) = meds.T1SCcmidopany(find(meds.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing med data \n'))
            m(sub,1) = 0;
        end
        % symptoms
        % if ~isempty(find(trilevel.PID==id{sub,1}))
        %     s(sub,1) = trilevel.T2gendis(find(trilevel.PID==id{sub,1}));
        %     s(sub,2) = trilevel.T2anhedon(find(trilevel.PID==id{sub,1}));
        %     s(sub,3) = trilevel.T2fears(find(trilevel.PID==id{sub,1}));
        % else
        %     fprintf(strcat(num2str(id{sub,1}),': missing clinical data \n'))
        %     s(sub,1) = NaN;
        %     s(sub,2) = NaN;
        %     s(sub,3) = NaN;
        % end
        % symptoms longitudinal
        if ~isempty(find(gd_long.PID==id{sub,1}))
            slong(sub,1) = gd_long.slope(find(gd_long.PID==id{sub,1}));
            slong(sub,2) = anh_long.slope(find(anh_long.PID==id{sub,1}));
            slong(sub,3) = fear_long.slope(find(fear_long.PID==id{sub,1}));
            
            sint(sub,1) = gd_long.intercept(find(gd_long.PID==id{sub,1}));
            sint(sub,2) = anh_long.intercept(find(anh_long.PID==id{sub,1}));
            sint(sub,3) = fear_long.intercept(find(fear_long.PID==id{sub,1}));
        else
            fprintf(strcat(num2str(id{sub,1}),': missing clinical data \n'))
            slong(sub,1) = NaN;
            slong(sub,2) = NaN;
            slong(sub,3) = NaN;
            sint(sub,1) = NaN;
            sint(sub,2) = NaN;
            sint(sub,3) = NaN;
        end
        
    end
    
    % site, age, sex, racial and ethnic identity, and participant reports
    % of annual family income - covariates from greg's paper.
    cti_temp = table2array([cti_final(:,6:9),cti_final(:,14:17),cti_final(:,22:25),cti_final(:,30:33),cti_final(:,38:41),cti_final(:,46:49)]);
    cti_temp = sum(cti_temp,2); cti_temp(cti_temp>0)=1;
    cti_traum = cti_temp;
    site = cell2mat(id); site(site<20000) = 0; site(site>20000) = 1; site = array2table(site);
    slong = array2table(slong); slong.Properties.VariableNames = {'longGeneralDistress','longAnhedoniaApprehension','longFears'};
    s = array2table(sint); s.Properties.VariableNames = {'GeneralDistress','AnhedoniaApprehension','Fears'};
    m = array2table(m); m.Properties.VariableNames = {'medication'};
    i = array2table(i); i.Properties.VariableNames = {'inflammation'};
    %d = array2table(d); d.Properties.VariableNames = dem.Properties.VariableNames;
    
    final_outcomes = [slong,s,m,i,d,site,array2table(traum1),array2table(traum2),array2table(inc),array2table(cti_traum)];
    
    % pull together the variables that you're actually going to use
    regressors = [final_outcomes.inflammation,final_outcomes.sex,final_outcomes.site,final_outcomes.medication,...
        final_outcomes.race,final_outcomes.ethnicity,final_outcomes.inc,final_outcomes.cti_traum,final_outcomes.GeneralDistress,final_outcomes.AnhedoniaApprehension,...
        final_outcomes.Fears,final_outcomes.longGeneralDistress,final_outcomes.longAnhedoniaApprehension,final_outcomes.longFears];
    
    regressors = array2table(regressors);
    regressors.Properties.VariableNames = {'inflammation','sex','site','meds','race','ethnicity','inc','cti','GeneralDistress',...
        'Anhedonia','Fears','longGeneralDistress','longAnhedonia','longFears'};
end

if whole_brain_networks == 1
    %% whole brain relationships between brain, inflammation and brain, symptoms
    remove_missing_data = 1;
    % This will also be my sanity check about what networks are associated
    % with each of the seeds I've chosen to study. I need to output average
    % images for the subs still in the study (stored in fnames) and then
    % run simple regressions after converting r to z using atanh
    

    % you can change which seed you're looking at here
    if overwrite_nii == 1
        if mid == 1
            for f = 1:length(fnames)
                load(fnames{f});
                pid=fnames{f}(1:5);
                final_brain=fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/beta_series/sub-10001/ses-2/gain_contrasts/run-1/beta_0001.nii');
                final_brain.dat = atanh(final_corr_vs_whole); % vs amyg acc ofc
                final_brain.fullpath = strcat(pid,'_vs_antgain.nii');
                write(final_brain,'overwrite')
            end
        end
        if rest == 1
            for f = 1:length(fnames)
                load(fnames{f});
                pid=fnames{f}(1:5);
                final_brain=fmri_data('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/beta_series/sub-10001/ses-2/gain_contrasts/run-1/beta_0001.nii');
                final_brain.dat = atanh(corr_vs_wholebrain); % vs amyg acc ofc
                final_brain.fullpath = strcat(pid,'_vs_rest.nii');
                write(final_brain,'overwrite')
            end
        end
    end
    % the sections commented out below were run once and then I saved the
    % final data object so I wouldn't have to regenerate it. It's still
    % there so I can go back and recreate if needed.
    full_fnames = filenames(fullfile(basedir,strcat('/*',region_name_for_wholebrain_analysis,'*nii')));
    final_brain = fmri_data(full_fnames);
    
    final_brain.X = ones(length(full_fnames),1); 
    %final_brain.X = [regressors.inflammation.*regressors.cti,regressors.cti,regressors.inflammation,regressors.sex,regressors.site,regressors.meds,regressors.race,regressors.ethnicity,regressors.inc];
    if remove_missing_data == 1
        % remove NaNs related to inflammation
        final_brain.dat(:, isnan(regressors.inflammation(:,1))) = [];
        full_fnames(isnan(regressors.inflammation(:,1)),:)=[];
        final_brain.X(isnan(regressors.inflammation(:,1)),:) = [];
        regressors(isnan(regressors.inflammation(:,1)),:) = [];
        
        % remove NaNs related to family income
        final_brain.dat(:, isnan(regressors.inc(:,1))) = [];
        full_fnames(isnan(regressors.inc(:,1)),:)=[];
        final_brain.X(isnan(regressors.inc(:,1)),:) = [];
        regressors(isnan(regressors.inc(:,1)),:) = [];
                   
        % remove NaNs related to symptoms
        final_brain.dat(:, isnan(regressors.longGeneralDistress(:,1))) = [];
        full_fnames(isnan(regressors.longGeneralDistress(:,1)),:)=[];
        final_brain.X(isnan(regressors.longGeneralDistress(:,1)),:) = [];
        regressors(isnan(regressors.longGeneralDistress(:,1)),:) = [];

        % remove NaNs related to cti
        final_brain.dat(:, isnan(regressors.cti(:,1))) = [];
        full_fnames(isnan(regressors.cti(:,1)),:)=[];
        final_brain.X(isnan(regressors.cti(:,1)),:) = [];
        regressors(isnan(regressors.cti(:,1)),:) = [];
    end
    statobj = regress(final_brain);
    threshobj = threshold(statobj.t,0.001,'unc','k',5);
    orthviews(select_one_image(threshobj,1))
    table(select_one_image(threshobj,1))
%     thresh_obj = final_brain;
%     thresh_obj.dat = mean(thresh_obj.dat,2);
%     thresh_obj.dat(thresh_obj.dat<0.2)=0;
%     orthviews(thresh_obj)
    
end    


%% Aim 1: seed to seed version of analyses above
if linear_seed_to_seed == 1
    remove_missing_data = 1;
    % the sections commented out below were run once and then I saved the
    % final data object so I wouldn't have to regenerate it. It's still
    % there so I can go back and recreate if needed.
    s2sfnames = filenames(fullfile(basedir,strcat('/*',region_name_for_wholebrain_analysis,'*nii')));
    
    if mid == 1
        s2sfnames(185) = [];
    end


    final_brain = fmri_data(s2sfnames);
    if exist('seed_to_seed_results','dir') == 0
        mkdir seed_to_seed_results/
    end
    
    cd('seed_to_seed_results/')
    amyg = fmri_data(filenames(fullfile(seeddir,'*Amygdala*')));
    vs = fmri_data(filenames(fullfile(seeddir,'VS*Rew.nii')));
    acc = fmri_data(filenames(fullfile(seeddir,'*ACC.nii')));
    ofc = fmri_data(filenames(fullfile(seeddir,'OFC*.nii')));
    % apply masks to data
    amyg_data = extract_roi_averages(final_brain,amyg);
    vs_data = extract_roi_averages(final_brain,vs);
    ofc_data = extract_roi_averages(final_brain,ofc);
    acc_data = extract_roi_averages(final_brain,acc);

    seed2seed = array2table([amyg_data.dat,vs_data.dat,ofc_data.dat,acc_data.dat]); 
    seed2seed.Properties.VariableNames = {'amygdala','vs','ofc','acc'};

     % combine with regressors
    alldata = [regressors,seed2seed];

    if remove_missing_data == 1
        % remove NaNs related to inflammation
        s2sfnames(isnan(alldata.inflammation(:,1)),:)=[];
        alldata(isnan(alldata.inflammation(:,1)),:) = [];
        
        % remove NaNs related to family income
        s2sfnames(isnan(alldata.inc(:,1)),:)=[];
        alldata(isnan(alldata.inc(:,1)),:) = [];
            
%         % remove NaNs related to symptoms
%         s2sfnames(isnan(alldata.longGeneralDistress(:,1)),:)=[];
%         alldata(isnan(alldata.longGeneralDistress(:,1)),:) = [];
% 
%         % remove NaNs related to cti
%         s2sfnames(isnan(alldata.cti(:,1)),:)=[];
%         alldata(isnan(alldata.cti(:,1)),:) = [];
    end
    
   

    if strcmp(region_name_for_wholebrain_analysis,'acc')
        mdlamyg = fitlm(alldata,'amygdala  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlofc = fitlm(alldata,'ofc  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlvs = fitlm(alldata,'vs  ~ inflammation + site + sex + meds + race + inc + ethnicity')
    elseif strcmp(region_name_for_wholebrain_analysis,'ofc')
        mdlamyg = fitlm(alldata,'amygdala  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlacc = fitlm(alldata,'acc  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlvs = fitlm(alldata,'vs  ~ inflammation + site + sex + meds + race + inc + ethnicity')
    elseif strcmp(region_name_for_wholebrain_analysis,'vs')
        mdlamyg = fitlm(alldata,'amygdala  ~ inflammation + site + sex + race + inc + ethnicity')
        mdlacc = fitlm(alldata,'acc  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlofc = fitlm(alldata,'ofc  ~ inflammation + site + sex + meds + race + inc + ethnicity')
    elseif strcmp(region_name_for_wholebrain_analysis,'amyg')
        mdlvs = fitlm(alldata,'vs  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
        mdlacc = fitlm(alldata,'acc  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
        mdlofc = fitlm(alldata,'ofc  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
    end

    if mediation_seed_to_seed == 1
        % output doesn't save super nicely so I'm not going to bother. Will
        % run each line manually
        save_all_mdls = 0;
        % only writing code for VS connections to minimize tests
        if save_all_mdls == 1
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.amygdala,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoAmygdala'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.amygdala,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoAmygdala'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.amygdala,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoAmygdala'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.acc,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoACC'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.acc,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoACC'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.acc,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoACC'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.ofc,'names',{'X:inflammation' 'Y:GeneralDistress' 'M:VStoOFC'},'verbose','dosave','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.ofc,'names',{'X:inflammation' 'Y:Anhedonia' 'M:VStoOFC'},'verbose','dosave','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.ofc,'names',{'X:inflammation' 'Y:Fears' 'M:VStoOFC'},'verbose','dosave','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
        else
            keyboard
        end
    end

    % AIM 3 MODERATED MEDIATION
    if moderated_mediation_seed_to_seed == 1
        % output doesn't save super nicely so I'm not going to bother. Will
        % run each line manually
        save_all_mdls = 0;
        % only writing code for VS connections to minimize tests
        if save_all_mdls == 1
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.amygdala.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoAmygdala*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.amygdala.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoAmygdala*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.amygdala.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoAmygdala*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.acc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoACC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.acc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoACC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.acc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:VStoACC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.ofc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:GeneralDistress' 'M:VStoOFC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.ofc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Anhedonia' 'M:VStoOFC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.ofc.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Fears' 'M:VStoOFC*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
        else
            keyboard
        end
    end
end

%% whole brain mediation: aim 2
% here I extend the roi analyses presented in aims 1 and 2 to understand if
% there is anywhere in the brain that appears to systematically covary with
% variables of interest in the context of mediation

if whole_brain_mediation_analysis == 1
    remove_missing_data = 1;
    full_fnames = filenames(fullfile(basedir,strcat('/*',region_name_for_wholebrain_analysis,'*nii')));

    if mid == 1
        full_fnames(185) = [];
    end
    
    if remove_missing_data == 1
        % remove NaNs related to inflammation
        full_fnames(isnan(regressors.inflammation(:,1)),:)=[];
        regressors(isnan(regressors.inflammation(:,1)),:) = [];
        
        % remove NaNs related to family income
        full_fnames(isnan(regressors.inc(:,1)),:)=[];
        regressors(isnan(regressors.inc(:,1)),:) = [];
            
        % remove NaNs related to symptoms
        full_fnames(isnan(regressors.longGeneralDistress(:,1)),:)=[];
        regressors(isnan(regressors.longGeneralDistress(:,1)),:) = [];

        % remove NaNs related to cti
        full_fnames(isnan(regressors.cti(:,1)),:)=[];
        regressors(isnan(regressors.cti(:,1)),:) = [];
    end

    final_brain = fmri_data(full_fnames);
    % set up directories
    if exist("mediation_results/",'dir') == 0
        mkdir mediation_results
    end

    cd mediation_results

    if exist('general_distress','dir') == 0   
        mkdir anhedonia; mkdir fears; mkdir general_distress
    end
        
    names = {'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoWholebrain'};
    mask = which('gray_matter_mask.img');
    
    cd(fullfile(basedir,'/mediation_results/general_distress'))

    resultsgendis = mediation_brain(regressors.inflammation, regressors.longGeneralDistress, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);    
    
    cd(fullfile(basedir,'/mediation_results/anhedonia'))

    resultsanhed = mediation_brain(regressors.inflammation, regressors.longAnhedonia, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);

    cd(fullfile(basedir,'/mediation_results/fears'))

    resultsfears = mediation_brain(regressors.inflammation, regressors.longFears, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);

    % use this code to change around folders and check for results
    mediation_brain_results('ab', 'thresh', ...
    0.001, ... % can insert for cluster thresh: , 'size', 5
    'whitebackground','slices', 'tables', 'names', 'save');
end


% AIM 3 MODERATED MEDIATION
if whole_brain_moderated_mediation == 1
    remove_missing_data = 1;
    full_fnames = filenames(fullfile(basedir,strcat('/*',region_name_for_wholebrain_analysis,'*nii')));

    if mid == 1
        full_fnames(185) = [];
    end
    
    if remove_missing_data == 1
        % remove NaNs related to inflammation
        full_fnames(isnan(regressors.inflammation(:,1)),:)=[];
        regressors(isnan(regressors.inflammation(:,1)),:) = [];
        
        % remove NaNs related to family income
        full_fnames(isnan(regressors.inc(:,1)),:)=[];
        regressors(isnan(regressors.inc(:,1)),:) = [];
            
        % remove NaNs related to symptoms
        full_fnames(isnan(regressors.longGeneralDistress(:,1)),:)=[];
        regressors(isnan(regressors.longGeneralDistress(:,1)),:) = [];

        % remove NaNs related to cti
        full_fnames(isnan(regressors.cti(:,1)),:)=[];
        regressors(isnan(regressors.cti(:,1)),:) = [];
    end

    final_brain = fmri_data(full_fnames);
    % set up directories
    if exist("moderated_mediation_results/",'dir') == 0
        mkdir moderated_mediation_results
    end

    cd moderated_mediation_results

    if exist('general_distress','dir') == 0   
        mkdir anhedonia; mkdir fears; mkdir general_distress
    end
        
    names = {'X:inflammation' 'Y:longitudinalsymptom' 'M:VStoWholebrain'};
    mask = which('gray_matter_mask.img');
    
    cd(fullfile(basedir,'/moderated_mediation_results/general_distress'))

    resultsgendis = mediation_brain(regressors.inflammation, regressors.longGeneralDistress, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);    
    
    cd(fullfile(basedir,'/moderated_mediation_results/anhedonia'))

    resultsanhed = mediation_brain(regressors.inflammation, regressors.longAnhedonia, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);

    cd(fullfile(basedir,'/moderated_mediation_results/fears'))

    resultsfears = mediation_brain(regressors.inflammation, regressors.longFears, final_brain.fullpath,'names',names,'mask', mask,'noverbose','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);

    % use this code to change around folders and check for results
    mediation_brain_results('ab', 'thresh', ...
    0.001, ... % can insert for cluster thresh: , 'size', 5
    'whitebackground','slices', 'tables', 'names', 'save');

end

% no significant clusters are found for the indirect effect from
% inflammation through to longitudinal symptoms. I could imagine doing a
% moderated mediation here as well. Though given the lack of findings in
% the seed to seed version of analyses, that doesn't feel like a smart
% move.
    
%% hyperalignment 

if hyper == 1
     
    % only one of these can be turned 'on' at a time. Composites should be
    % used when you're 1) trying to optimize speed or 2) want to get voxelwise
    % translation. Basically, this feeds hyperalignment matrices that yield a
    % 162x163 rotation matrix, the less computationally expensive route. On the
    % other hand, if you need voxelwise rotation information, then grab
    % rotation. It will take longer and it will be important to change the perm
    % loop to be something like 1:10. 1:1000 will crash your computer
    composites = 1;
    rotation = 0;
    % to optimize signal for plotting, I'm going to feed hyperalignment
    % controls first and psychosis second. I'll randomize the order of subjects
    % across 10 permutations, but am trying to stack the deck to build a
    % hyperaligned signature of psychosis in this sample
    pure_random = 1;
    % when rotation = 1, this value should be small (like 10). The rotation
    % method mentioned above is super expensive computationally. If you're
    % going for composites then feel free to set this up at 1000. That
    % shouldn't crash your computer and you'll want the extra iterations. This
    % is especially true if you're doing a pure random test, which will shuffle
    % subjects in all kinds of different orders.
    permutation_length = 100;
    
    for f = 1:length(fnames)
        load(fnames{f});
        % this will load three variables. corr_vs_to_wholebrain, hypmat,
        % temp_brain. the one you want for this set of analyses is hypmat
        % which will be a VS voxel x 300 ROI masked connectivity matrix.
        % Available variables are: final_corr_acc_voxel_to_region
        % final_corr_ofc_voxel_to_region final_corr_amyg_voxel_to_region
        % final_corr_vs_voxel_to_region 
        if mid == 1
            unaligned_mats{f} = final_corr_vs_voxel_to_region;
        end
        if rest == 1
            unaligned_mats{f} = corr_vs_voxel_to_region;
        end
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
        end
        fprintf(strcat('Im on permutation: ',num2str(perm),'\n'))
        curr_filename1 = fullfile(basedir,'/hyperalignment_results/',strcat('amyg_data_perm_',num2str(perm),'.mat'));
        curr_filename2 = fullfile(basedir,'/hyperalignment_results/',strcat('amyg_transforms_perm_',num2str(perm),'.mat'));
        save(curr_filename1, 'aligned_data_resorted', 'unaligned_mats');
        save(curr_filename2,'transforms_resorted', '-v7.3')
        % this bit calculates and stores composite values       

    end
%     cd(fullfile(basedir,'hyperalignment_results'))
%     if overwrite_hyper == 1
%         save predict_hyperalignment_composite_perm.mat aligned_data_resorted unaligned_mats transforms_resorted -v7.3  
%     end
end

    
if hyper_analyze == 1
    cd(fullfile(basedir,'/hyperalignment_results'))
    
    % need to combine across 100 permutations
    analyze_fnames = filenames(fullfile(strcat(region_name_for_wholebrain_analysis,'_data_perm*')));
    
    % this takes a long time. I'm going to save each roi's version of this
    % output so I can simply load it in when I want to relate new variables
    % to previously processed brain data. I'm going to save the output in a
    % shortcuts/ folder in hyperalignment results
    if remake_multivar_input == 1
        for i = 1:100
            load(analyze_fnames{i})
            for sub = 1:length(aligned_data_resorted)
                if i == 1
                    final_aligned_data.(strcat('perm',num2str(i)))(:,sub) = aligned_data_resorted{sub}(:);
                    final_unaligned_data(:,sub) = unaligned_mats{sub}(:);
                else                
                    final_aligned_data(:,sub).(strcat('perm',num2str(i))) = aligned_data_resorted{sub}(:);
                end
            end
        end
    else
        fname = filenames(fullfile('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/hyperalignment_results/shortcuts/',strcat(region_name_for_wholebrain_analysis,'*')));
        load(fname{1});
    end
    keyboard
    if do_pls_regress == 1

        % easier to define an outcome here which will then be plugged in
        % throughout the document
        outcome = regressors.inflammation(:,1);

        % remove NaNs related to outcome
        final_unaligned_data(:, isnan(outcome(:,1))) = [];
        for i = 1:100
            final_aligned_data.(strcat('perm',num2str(i)))(:, isnan(outcome(:,1))) = [];
        end
        outcome(isnan(outcome(:,1))) = [];    


        % saving this code in case I want to pivot back to plsregression
        [uXL,uYL,uXS,uYS,uBETA,uPCTVAR,uMSE,~] = plsregress(final_unaligned_data',outcome,10,'CV',10);

        for i=1:length(aligned_data_resorted)
            aligned_x = final_aligned_data.(strcat('perm',num2str(i)))(:,:);
            fprintf(strcat('Im on permutation: ',num2str(perm),'\n'))
            [aXL{perm},aYL{perm},aXS{perm},aYS{perm},aBETA{perm},aPCTVAR{perm},aMSE{perm},~] = plsregress(aligned_x',outcome,10,'CV',10);    
            clear aligned_x unaligned_x            
        end

        for i = 1:length(aMSE)
            aMSEall(i,:) = aMSE{i}(2,:);
            aPCTVARall(i,:) = aPCTVAR{i}(2,:);
        end
        finalMSEa = mean(mean(aMSEall),2);
        finalMSEu = mean(uMSE(2,:));
        avgPVEa = mean(aPCTVARall);
        
        load in atlas to highlight regions connected with VS in each of the
        above components
        

        figure(); plot(1:length(uYL),cumsum(100*avgPVEa),'-bo');
        title('Aligned data')
        xlabel('Number of PLS components');
        ylabel('Percent Variance Explained in inflammation');

        figure(); plot(1:length(uYL),cumsum(100*uPCTVAR(2,:)),'-bo');
        title('Unaligned data')
        xlabel('Number of PLS components');
        ylabel('Percent Variance Explained in inflammation');
    end
    for perm = 1:length(transforms_resorted)
        for sub = 1:length(transforms_resorted{1})
            % magnitude of translation (distance metric)
            
            all_translation(sub,:,perm) = (sum((transforms_resorted{perm}{sub}.c.^2),2)).^0.5;      
            % magnitude of rotation (distance metric)
            
            all_rotation(sub,:,perm) = (sum(abs(transforms_resorted{perm}{sub}.T),2));
            % magnitude of scaling (single value)
            all_scale(sub,perm) = transforms_resorted{perm}{sub}.b;
            
            % determinant of rotation matrix (were data rotated with
            % reflection (det(T) = -1, or not?
            all_reflection(sub,perm) = det(transforms_resorted{perm}{sub}.T);          
        end
    end
    % analyze transformation matrices. first total translation and rotation
    sub_total_translation = mean(zscore(mean(all_translation,3)),2);
    sub_total_rotation = mean(zscore(mean(all_rotation,3)),2);
    sub_total_scale = zscore(mean(all_scale,2));
    hyper_mdl_inputs = [sub_total_translation,sub_total_rotation,sub_total_scale];
    hyper_mdl_inputs = array2table(hyper_mdl_inputs); hyper_mdl_inputs.Properties.VariableNames = {'translation','rotation','scale'};
    mdl_input = [regressors,hyper_mdl_inputs];

    mdl_inf_trans = fitlm(mdl_input,'inflammation ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_inf_rot = fitlm(mdl_input,'inflammation ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_inf_scale = fitlm(mdl_input,'inflammation ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_gd_trans = fitlm(mdl_input,'GeneralDistress ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_gd_rot = fitlm(mdl_input,'GeneralDistress ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_gd_scale = fitlm(mdl_input,'GeneralDistress ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_anh_trans = fitlm(mdl_input,'Anhedonia ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_anh_rot = fitlm(mdl_input,'Anhedonia ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_anh_scale = fitlm(mdl_input,'Anhedonia ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_fear_trans = fitlm(mdl_input,'Fears ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_fear_rot = fitlm(mdl_input,'Fears ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_fear_scale = fitlm(mdl_input,'Fears ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_lgd_trans = fitlm(mdl_input,'longGeneralDistress ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lgd_rot = fitlm(mdl_input,'longGeneralDistress ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_lgd_scale = fitlm(mdl_input,'longGeneralDistress ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_lanh_trans = fitlm(mdl_input,'longAnhedonia ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lanh_rot = fitlm(mdl_input,'longAnhedonia ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_lanh_scale = fitlm(mdl_input,'longAnhedonia ~ scale + race + ethnicity + sex + meds + inc + site');

    mdl_lfear_trans = fitlm(mdl_input,'longFears ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lfear_rot = fitlm(mdl_input,'longFears ~ rotation + race + ethnicity + sex + meds + inc + site');
    mdl_lfear_scale = fitlm(mdl_input,'longFears ~ scale + race + ethnicity + sex + meds + inc + site');

end

if compare_aligned_pls_loadings == 1
    cd('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices/anticipation_combined_runs_with_hyp/hyperalignment_results')
    load('predict_hyperalignment_composite_perm.mat')
    load('group_level_regressors.mat')

    % the goal of this section will be to generate some quick similar
    % matrices of the feature loadings with respect to inflammation and all
    % 3 symptoms. The goal is to identify features (VS-wholebrain
    % connections) that appear relatively specific to each 
    inflam = load('predict_inflammation_hyperalignment_model_results.mat');
    gd = load('predict_gd_hyperalignment_model_results.mat');
    anhed = load('predict_anh_hyperalignment_model_results.mat');
    fears = load('predict_fears_hyperalignment_model_results.mat');
    for sub = 1:134
        for perm = 1:100
            if perm < 2
                loadings_inflam(:,:) = reshape(inflam.aXL{perm}(:,1),[300,514]);
                loadings_gd(:,:) = reshape(gd.aXL{perm}(:,1),[300,514]);
                loadings_anhed(:,:) = reshape(anhed.aXL{perm}(:,1),[300,514]);
                loadings_fears(:,:) = reshape(fears.aXL{perm}(:,1),[300,514]);
                % rotate loadings back into brain space
                
                rotated_c1_inflam(:,:,sub) = transforms_resorted{perm}{sub}.T(:,:)' * loadings_inflam(:,:);
                rotated_c1_gd(:,:,sub) = transforms_resorted{perm}{sub}.T(:,:)' * loadings_gd(:,:);
                rotated_c1_anhed(:,:,sub) = transforms_resorted{perm}{sub}.T(:,:)' * loadings_anhed(:,:);
                rotated_c1_fears(:,:,sub) = transforms_resorted{perm}{sub}.T(:,:)' * loadings_fears(:,:);
            else
                loadings_inflam(:,:) = reshape(inflam.aXL{perm}(:,1),[300,514]);
                loadings_gd(:,:) = reshape(gd.aXL{perm}(:,1),[300,514]);
                loadings_anhed(:,:) = reshape(anhed.aXL{perm}(:,1),[300,514]);
                loadings_fears(:,:) = reshape(fears.aXL{perm}(:,1),[300,514]);
                % rotate loadings back into brain space
                rotated_c1_inflam(:,:,sub) = (rotated_c1_inflam(:,:,sub) + transforms_resorted{perm}{sub}.T(:,:)' * loadings_inflam(:,:))./2;
                rotated_c1_gd(:,:,sub) = (rotated_c1_gd(:,:,sub) + transforms_resorted{perm}{sub}.T(:,:)' * loadings_gd(:,:))./2;
                rotated_c1_anhed(:,:,sub) = (rotated_c1_anhed(:,:,sub) + transforms_resorted{perm}{sub}.T(:,:)' * loadings_anhed(:,:))./2;
                rotated_c1_fears(:,:,sub) = (rotated_c1_fears(:,:,sub) + transforms_resorted{perm}{sub}.T(:,:)' * loadings_fears(:,:))./2;
            end
        end
    end
        

    total_rotated_c1_inflam = zscore(mean(sum(sqrt(sum(rotated_c1_inflam.^2,3)),2),3));
    total_rotated_c1_gd = zscore(mean(sum(sqrt(sum(rotated_c1_gd.^2,3)),2),3));
    total_rotated_c1_anhed = zscore(mean(sum(sqrt(sum(rotated_c1_anhed.^2,3)),2),3));
    total_rotated_c1_fears = zscore(mean(sum(sqrt(sum(rotated_c1_fears.^2,3)),2),3)); 

    atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');

    brain_connections_inflam = atl;
    brain_connections_gd = atl;
    brain_connections_anhed = atl;
    brain_connections_fears = atl;

    for r = 1:length(total_rotated_c1_inflam)
        brain_connections_inflam.dat(atl.dat==r) = total_rotated_c1_inflam(r);
        brain_connections_gd.dat(atl.dat==r) = total_rotated_c1_gd(r);
        brain_connections_anhed.dat(atl.dat==r) = total_rotated_c1_anhed(r);
        brain_connections_fears.dat(atl.dat==r) = total_rotated_c1_fears(r);
    end

    % pseudo thresholded brains: top 2.5% of connections

    brain_connections_inflam.dat(brain_connections_inflam.dat<2) = 0;
    brain_connections_gd.dat(brain_connections_gd.dat<2) = 0;
    brain_connections_anhed.dat(brain_connections_anhed.dat<2) = 0;
    brain_connections_fears.dat(brain_connections_fears.dat<2) = 0;
    
    r_brain_connections_inflam = region(brain_connections_inflam);
    r_brain_connections_gd = region(brain_connections_gd);
    r_brain_connections_anhed = region(brain_connections_anhed);
    r_brain_connections_fears = region(brain_connections_fears); 
    

end

if look_at_multivariate_rotation_estimates == 1
    cd('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/beta_series/final_matrices/anticipation_combined_runs_with_hyp/hyperalignment_results')
    load('predict_hyperalignment_composite_perm.mat')
    load('group_level_regressors.mat')
    atl = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/300ROIatlas/300_ROI_Set/ROIs_300inVol_MNI.nii');

    for perm = 1:length(transforms_resorted)
        if perm == 1
            for sub = 1:length(transforms_resorted{1})
                final_rots(sub,:) = sum(abs(transforms_resorted{perm}{sub}.T(:,:)),2);
                final_trans(sub,:) = sum(abs(transforms_resorted{perm}{sub}.c(:,:)),1);
            end
        else
            for sub = 1:length(transforms_resorted{perm})
                final_rots(sub,:) = sum(abs(transforms_resorted{perm}{sub}.T(:,:)),2);
                final_trans(sub,:) = sum(abs(transforms_resorted{perm}{sub}.c(:,:)),1)';
            end
        end
    end
    
    final_rots = zscore(final_rots);
    final_trans = zscore(final_trans);

    % series of models to test for differences in translation or rotation
    % that relate to inflammation

    for r = 1:size(final_rots,2)
        rot_temp = fitlm([regressors.inflammation,regressors.sex,regressors.site,regressors.inc, regressors.meds,...
            regressors.race,regressors.ethnicity],final_rots(:,r));
        tran_temp = fitlm([regressors.inflammation,regressors.sex,regressors.site,regressors.inc, regressors.meds,...
            regressors.race,regressors.ethnicity],final_trans(:,r));
        rot_beta(r,1) = rot_temp.Coefficients.tStat(2);
        tran_beta(r,1) = tran_temp.Coefficients.tStat(2);
        rot_p(r,1) = rot_temp.Coefficients.pValue(2);
        tran_p(r,1) = tran_temp.Coefficients.pValue(2);
    end

    % fdr correct images using Benjamini & Hochberg 2005
    [~,tidx] = sort(tran_p(:,1),'ascend');
    sorted_translation = tran_p(tidx,1);
    [~,t2idx]=sort(tidx);

    [~,ridx] = sort(rot_p(:,1),'ascend');
    sorted_rotation = rot_p(ridx,:);
    [~,r2idx]=sort(ridx);
    

    % FDR<0.1 correction
    for rank = 1:length(tran_p)
       fdr_p(rank,1) = (rank/300) * (0.1);
    end
 
    sorted_translation = [sorted_translation,fdr_p];
    sorted_rotation = [sorted_rotation,fdr_p];

    translation_maps_final = sorted_translation(t2idx,:);
    rotation_maps_final = sorted_rotation(r2idx,:);
    
    % changing this to be a 0.001 uncorrected threshold
    sig_tran = find(translation_maps_final(:,1)<0.001);%translation_maps_final(:,2));
    sig_rot = find(rotation_maps_final(:,1)<0.001);%rotation_maps_final(:,2));

    sig_rot_brain = atl; sig_rot_brain.dat = zeros(length(sig_rot_brain.dat),1);
    sig_tran_brain = atl; sig_tran_brain.dat = zeros(length(sig_tran_brain.dat),1);

    for r = 1:length(sig_rot)
        sig_rot_brain.dat(atl.dat==sig_rot(r)) = rot_beta(sig_rot(r));
    end
    for r = 1:length(sig_tran)
        sig_tran_brain.dat(atl.dat==sig_tran(r)) = rot_beta(sig_tran(r));
    end
    % inflammation: reduced rotation insula (60) and heightened rotation dlpfc (124)
    % long general distress: nothing
    % long anhedonia: heightened translation frontal inferior triangularis (146),
    % heightened translation middle occipital (218)
    % long fears: nothing
    % intercept general distress:nothing
    % intercept anhedonia: reduced inferior parietal (220)

end


if mediate_with_pls_components == 1
    load('predict_inflammation_hyperalignment_model_results.mat')
    load('group_level_regressors.mat')
    % you should do mediation with one and two levels of mediators. This
    % will include the first two components extracted that appear to
    % account for quite a bit more of the variance in in inflammation
    for perm = 1:length(aXL)
        % generate average loadings
        m1(:,perm) = aXS{perm}(:,1);
        m2(:,perm) = aXS{perm}(:,2);
        m3(:,perm) = aXS{perm}(:,3);
    end
    m1 = mean(m1,2);
    m2 = mean(m2,2);
    m3 = mean(m3,2);

    m1u = uXS(:,1);
    m2u = uXS(:,2);
    m3u = uXS(:,3);
    
    
    % use factor scores in mediation model
    [paths, stats2] = mediation(regressors.inflammation, regressors.longAnhedonia, m1, 'plots', 'verbose', 'names', {'inflammation' 'long. Anhedonia' 'VS component 1'},'covs',[regressors.meds,regressors.sex,regressors.site,regressors.race,regressors.ethnicity,regressors.inc]);
   

    additional_reg = [m1,m2,m3]; additional_reg = array2table(additional_reg);
    additional_reg.Properties.VariableNames = {'Component1','Component2','Component3'};

    regressors = [regressors,additional_reg];
    fitlm(regressors,'Anhedonia ~ Component1 + sex + site + meds + race + ethnicity + inc') % significant relationship with Anhedonia
    fitlm(regressors,'GeneralDistress ~ Component1 + sex + site + meds + race + ethnicity + inc')
    fitlm(regressors,'Fears ~ Component1 + sex + site + meds + race + ethnicity + inc')

    fitlm(regressors,'Anhedonia ~ Component2 + sex + site + meds + race + ethnicity + inc') 
    fitlm(regressors,'GeneralDistress ~ Component2 + sex + site + meds + race + ethnicity + inc')
    fitlm(regressors,'Fears ~ Component2 + sex + site + meds + race + ethnicity + inc')

    fitlm(regressors,'Anhedonia ~ Component3 + sex + site + meds + race + ethnicity + inc') 
    fitlm(regressors,'GeneralDistress ~ Component3 + sex + site + meds + race + ethnicity + inc')
    fitlm(regressors,'Fears ~ Component3 + sex + site + meds + race + ethnicity + inc')
end
% to do
% brain mediation with factor loadings per voxel. use mediation_brain


% pull connections for unaligned data that contribute to relationship with
% each ourcome

if network_based_analyses == 1
    atl = fmri_data('ROIs_300inVol_MNI.nii');
    
    labels = readtable('ROIs_300inVol_MNI_allInfo.txt');
    
    region_names = readtable("ROIs_anatomicalLabels.txt"); % 0=cortexMid,1=cortexL,2=cortexR,3=hippocampus,4=amygdala,5=basalGanglia,6=thalamus,7=cerebellum
    
    cerebellum_nodes = flip(find(region_names.x0_cortexMid==7));
    for i = 1:length(cerebellum_nodes)
        labels(cerebellum_nodes(i),:) = [];
    end
    
    roi_ids_rew = find(contains(labels.netName,'Reward'));
    rewdat = atl; rewdat.dat(:,1) = 0;
    
    for i = 1:length(roi_ids_rew)
        rewdat.dat(atl.dat==roi_ids_rew(i)) = 1;
    end
    
    roi_ids_dm = find(contains(labels.netName,'Default'));
    
    dmdat = atl; dmdat.dat(:,1) = 0;
    
    for i = 1:length(roi_ids_dm)
        dmdat.dat(atl.dat==roi_ids_dm(i)) = 1;
    end
    
    roi_ids_fp = find(contains(labels.netName,'FrontoParietal'));
    
    fpdat = atl; fpdat.dat(:,1) = 0;
    
    for i = 1:length(roi_ids_fp)
        fpdat.dat(atl.dat==roi_ids_fp(i)) = 1;
    end


end
