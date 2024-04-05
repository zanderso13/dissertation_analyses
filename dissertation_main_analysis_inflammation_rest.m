% I think this one should always be equal to 1
prep_behavioral_data = 1;

% what task?
mid = 0; rest = 1;

region_name_for_wholebrain_analysis = 'vs'; % vs amyg acc ofc

whole_brain_networks = 0; overwrite_nii = 0; 

linear_seed_to_seed = 0; seeddir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/seeds';
mediation_seed_to_seed = 0; % linear seed to see must be on
moderated_mediation_seed_to_seed = 0;
network_based_analyses = 0; overwrite_networks = 0; mediation_networks = 0; 
moderated_mediation_networks = 0; % network based analyses must be on

whole_brain_mediation_analysis = 0; 
whole_brain_moderated_mediation = 0;

hyper= 0; % hardcoded which seed you're looking at
important_to_change = 1; % table of contents. This is the outcome. Need to remove NaN from the outcome you want to analyze

hyper_analyze = 1; do_pls_regress = 1; 

hyper_analyze_transforms = 1; % visualize these pls results on transformations using the visualize_pls_results option

visualize_pls_results = 0;

mediation_with_hyp_transform = 0;

if mid == 1
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats';
end

if rest == 1
    basedir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices';
end    

cd(basedir)

fnames = filenames('*.mat');
if mid==1   
    motion_exclude_temp = load(fullfile('/home/zaz3744/repo/dissertation_analyses/final_motion_exclusion.mat'));
    
    % add exclusions for falling asleep, problems with scanner task, etc
    admin_problems = {'10029','10319','10323','20025','20150','10230','10002','21057','21463'}; % incomplete scans, no timing files, asleep
    preprocessing_quality = {'10341','10247','10461','21025','20133','20309','20507','21046','21163','21675','20133','20309'}; % scan artificats, bad grayplots
    anatomical_quality = {'10141','20235','21417','21463'}; 
    exclude_winnings = {'10001','10094','10102','10125','10140','10423','10443','10461','10471','20032','20108','20564','20674','21111','21223'};
    temp_tot_list = [admin_problems,preprocessing_quality,anatomical_quality,exclude_winnings];
    
    for i = 1:length(temp_tot_list)
        motion_exclude_temp.pid_exclude_list(contains(motion_exclude_temp.pid_exclude_list(:,1),temp_tot_list{i}),:)=[];
    end
    
    pid_exclude_list = [temp_tot_list';motion_exclude_temp.pid_exclude_list(:,1)];
end
% if you want to do 0.3 fd cutoff for rest, you can actually just comment the final motion exclusion file out 
if rest==1   
    motion_exclude_temp = load(fullfile('/home/zaz3744/repo/dissertation_analyses/final_motion_exclusions_rest.mat'));
    % add exclusions for falling asleep, problems with scanner task, etc
    admin_problems = {'10004','21052','21163','21184'}; % incomplete scans, no timing files, asleep
    preprocessing_quality = {'10001','10002','10006','10008','10010','10034','10041','10059',...
        '10088','10090','10135','10272','10341','10422','20085','20123','20309','20464','21025',...
        '21463'}; % scan artificats, bad grayplots, field of view
    anatomical_quality = {'10141','20235','21417','21463'}; 
    temp_tot_list = [admin_problems,preprocessing_quality,anatomical_quality];
    
    for i = 1:length(temp_tot_list)
        motion_exclude_temp.pid_exclude_list(contains(motion_exclude_temp.pid_exclude_list(:,1),temp_tot_list{i}),:)=[];
    end
    
    pid_exclude_list = [temp_tot_list';motion_exclude_temp.pid_exclude_list(:,1)];
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
% if mid == 1
%     fnames(185) = [];
% end

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
            i(sub,1) = immune.T1BDicsavg(find(immune.PID==id{sub,1}));
            %i(sub,1) = immune.ZT1BDcrpLN(find(immune.PID==id{sub,1})); %
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
        %     s(sub,2) = trilevel.T2anhedon(find(trileveimportant_to_changel.PID==id{sub,1}));
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
    
    %final_brain.X = ones(length(full_fnames),1); 
    final_brain.X = [regressors.inflammation.*regressors.cti,regressors.cti,regressors.inflammation,regressors.sex,regressors.site,regressors.meds,regressors.race,regressors.ethnicity,regressors.inc];
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
    threshobj = threshold(statobj.t,0.001,'unc','k',10);
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

    final_brain = fmri_data(s2sfnames);
    if exist('seed_to_seed_results','dir') == 0
        mkdir seed_to_seed_results/
    end
    
    cd('seed_to_seed_results/')
    % old set of masks that come from too many different sources
%     amyg = fmri_data(filenames(fullfile(seeddir,'*Amygdala*')));
%     vs = fmri_data(filenames(fullfile(seeddir,'VS*Rew.nii')));
%     acc = fmri_data(filenames(fullfile(seeddir,'*ACC.nii')));
%     ofc = fmri_data(filenames(fullfile(seeddir,'OFC*.nii')));

    % generate new masks for each target region
    atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    
    % vs seitz idx = 246 247
    % amyg seitz idx = 244 245
    % ofc seitz idx = 105 111 116 117 118
    % acc seitz idx = 102 110 108 122 204 206 208

    vs = atl; vs.dat(vs.dat>0) = 0; vs.dat(atl.dat==244)=1; vs.dat(atl.dat==245)=1;

    amyg = atl; amyg.dat(amyg.dat>0) = 0; amyg.dat(atl.dat==244)=1; amyg.dat(atl.dat==245)=1;

    ofc = atl; ofc.dat(ofc.dat>0) = 0; ofc.dat(atl.dat==105)=1; ofc.dat(atl.dat==111)=1;
    ofc.dat(atl.dat==116)=1; ofc.dat(atl.dat==117)=1; ofc.dat(atl.dat==118)=1;

    acc = atl; acc.dat(acc.dat>0) = 0; acc.dat(atl.dat==102)=1; acc.dat(atl.dat==110)=1; ofc.dat(atl.dat==108)=1;
    acc.dat(atl.dat==122)=1; acc.dat(atl.dat==204)=1; acc.dat(atl.dat==206)=1; acc.dat(atl.dat==208)=1;

    % apply masks to data
    amyg_data = extract_roi_averages(final_brain,amyg);
    ofc_data = extract_roi_averages(final_brain,ofc);
    acc_data = extract_roi_averages(final_brain,acc);
    vs_data = extract_roi_averages(final_brain,vs);

    seed2seed = array2table([amyg_data.dat,vs_data.dat,ofc_data.dat,acc_data.dat]); 
    seed2seed.Properties.VariableNames = {'amygdala','vs','ofc','acc'};

     % combine with regressors
    alldata = [regressors,seed2seed];

    if remove_missing_data == 1
        % remove NaNs related to inflammation
        s2sfnames(isnan(alldata.inflammation(:,1)),:)=[];
        fnames(isnan(alldata.inflammation(:,1)),:) = [];
        alldata(isnan(alldata.inflammation(:,1)),:) = [];
        
        % remove NaNs related to family income
        s2sfnames(isnan(alldata.inc(:,1)),:)=[];
        fnames(isnan(alldata.inc(:,1)),:) = [];
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
        mdlamyg = fitlm(alldata,'amygdala  ~ inflammation*cti + site + sex + race + inc + ethnicity')
        mdlacc = fitlm(alldata,'acc  ~ inflammation*cti + site + sex + meds + race + inc + ethnicity')
        mdlofc = fitlm(alldata,'ofc  ~ inflammation*cti + site + sex + meds + race + inc + ethnicity')
    elseif strcmp(region_name_for_wholebrain_analysis,'amyg')
        mdlvs = fitlm(alldata,'vs  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
        mdlacc = fitlm(alldata,'acc  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
        mdlofc = fitlm(alldata,'ofc  ~ longGeneralDistress + site + sex + meds + race + inc + ethnicity')
    end
    %% AIM 2
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

    if network_based_analyses == 1
        atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
        
        labels = readtable('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI_allInfo.txt');
        
        region_names = readtable("/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_anatomicalLabels.txt"); % 0=cortexMid,1=cortexL,2=cortexR,3=hippocampus,4=amygdala,5=basalGanglia,6=thalamus,7=cerebellum
        
        cerebellum_nodes = flip(find(region_names.x0_cortexMid==7));
        for i = 1:length(cerebellum_nodes)
            labels(cerebellum_nodes(i),:) = [];
        end
        
        roi_ids_rew = find(contains(labels.netName,'Reward'));
        rewdat = atl; rewdat.dat(:,1) = 0;
        
        for i = 1:length(roi_ids_rew)
            rewdat.dat(atl.dat==roi_ids_rew(i)) = 1;
        end

        roi_ids_co = find(contains(labels.netName,'CinguloOpercular'));
        codat = atl; codat.dat(:,1) = 0;

        for i = 1:length(roi_ids_co)
            codat.dat(atl.dat==roi_ids_co(i)) = 1;
        end
        
        for i = 1:length(roi_ids_rew)
            rewdat.dat(atl.dat==roi_ids_rew(i)) = 1;
        end
        
        roi_ids_dm = find(contains(labels.netName,'Default'));
        
        dmdat = atl; dmdat.dat(:,1) = 0;
        dmdatall = atl; dmdatall.dat(:,1) = 0;
        
        for i = 1:length(roi_ids_dm)
            dmdat.dat(atl.dat==roi_ids_dm(i)) = 1;
            dmdatall.dat(atl.dat==roi_ids_dm(i)) = roi_ids_dm(i);
        end
        
        roi_ids_fp = find(contains(labels.netName,'FrontoParietal'));
        
        fpdat = atl; fpdat.dat(:,1) = 0;
        fpdatall = atl; fpdatall.dat(:,1) = 0;
        
        for i = 1:length(roi_ids_fp)
            fpdat.dat(atl.dat==roi_ids_fp(i)) = 1;
            fpdatall.dat(atl.dat==roi_ids_fp(i)) = roi_ids_fp(i);
        end

        roi_ids_s = find(contains(labels.netName,'Salience'));
        
        saldat = atl; saldat.dat(:,1) = 0;
        saldatall = atl; saldatall.dat(:,1) = 0;
        
        for i = 1:length(roi_ids_s)
            saldat.dat(atl.dat==roi_ids_s(i)) = 1;
            saldatall.dat(atl.dat==roi_ids_s(i)) = roi_ids_s(i);
        end

        if overwrite_networks == 1
            if rest == 1
                cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices')
            end
            if mid == 1
                cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats')
            end
            for i = 1:length(fnames)
                fprintf(strcat('Working on:',fnames{i}(1:5),'\n'))
                load(fnames{i});
                if i==1
                    if mid==1
                        final_seitz_mat = seitz_mat;
                    end
                    if rest==1
                        final_seitz_mat{i} = seitz_mat';
                    end
                else
                    if mid==1
                        final_seitz_mat(:,:,i) = seitz_mat;
                    end
                    if rest==1
                        final_seitz_mat{i}(:,:) = seitz_mat';
                    end
                end
            end
            cd('seed_to_seed_results')
            save seitzman_time_series.mat final_seitz_mat
            
        else
            if rest == 1
                load('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices/seed_to_seed_results/seitzman_time_series.mat')
            end

            if mid == 1
                load('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/seed_to_seed_results/seitzman_time_series.mat')
            end

        end
        if rest == 1
            for i = 1:length(final_seitz_mat) % cell array. so use length              
                xsubs_corr_mat(:,:,i) = corr(final_seitz_mat{i}');
            end
        end
        if mid == 1
            for i = 1:size(final_seitz_mat,3)              
                xsubs_corr_mat(:,:,i) = corr(final_seitz_mat(:,:,i)');
            end
        end
        
        %

        for i = 1:size(xsubs_corr_mat,3)
            temprew = xsubs_corr_mat(roi_ids_rew,roi_ids_rew,i);
            temprew(temprew==1)=[];
            connrew(i,1) = mean(temprew);

            tempco = xsubs_corr_mat(roi_ids_co,roi_ids_co,i);
            tempco(tempco==1)=[];
            connco(i,1) = mean(tempco);
            
            tempdmn = xsubs_corr_mat(roi_ids_dm,roi_ids_dm,i);
            tempdmn(tempdmn==1)=[];
            conndmn(i,1) = mean(tempdmn);


            tempfpn = xsubs_corr_mat(roi_ids_fp,roi_ids_fp,i);
            tempfpn(tempfpn==1)=[];
            connfpn(i,1) = mean(tempfpn);

            tempsal = xsubs_corr_mat(roi_ids_s,roi_ids_s,i);
            tempsal(tempsal==1)=[];
            connsal(i,1) = mean(tempsal);

            temprewtodmn = xsubs_corr_mat(roi_ids_rew,roi_ids_dm,i);
            temprewtodmn(temprewtodmn==1)=[];
            connrewtodmn(i,1) = mean(temprewtodmn(:));
            
        end
        networks_conn = [atanh(connrew),atanh(conndmn),atanh(connfpn),atanh(connsal),atanh(connrewtodmn),atanh(connco)];
        networks_conn = array2table(networks_conn);
        networks_conn.Properties.VariableNames = {'Reward','DefaultMode','FrontoParietal','Salience','RewardtoDMN','CinguloOpercular'};
        
        alldata = [alldata,networks_conn];

        mdlrew = fitlm(alldata,'Reward  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdldmn = fitlm(alldata,'DefaultMode  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlfpn = fitlm(alldata,'FrontoParietal  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlsal = fitlm(alldata,'Salience  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlrewdmn = fitlm(alldata,'RewardtoDMN  ~ inflammation + site + sex + meds + race + inc + ethnicity')
        mdlco = fitlm(alldata,'CinguloOpercular  ~ inflammation + site + sex + meds + race + inc + ethnicity')

        mdlmiddmn = fitlm(alldata,'GeneralDistress  ~ inflammation*DefaultMode + site + sex + meds + race + inc + ethnicity')
        mdlrestfpn = fitlm(alldata,'GeneralDistress  ~ inflammation*FrontoParietal + site + sex + meds + race + inc + ethnicity')

        mdlrew2 = fitlm(alldata,'Reward  ~ longGeneralDistress + longAnhedonia + longFears + site + sex + meds + race + inc + ethnicity')
        mdldmn2 = fitlm(alldata,'DefaultMode  ~ longGeneralDistress + longAnhedonia + longFears + site + sex + meds + race + inc + ethnicity')
        mdlfpn2 = fitlm(alldata,'FrontoParietal  ~ longGeneralDistress + longAnhedonia + longFears + site + sex + meds + race + inc + ethnicity')
        mdlsal = fitlm(alldata,'Salience  ~ longGeneralDistress + longAnhedonia + longFears + site + sex + meds + race + inc + ethnicity')
        mdlrewdmn2 = fitlm(alldata,'RewardtoDMN  ~ longGeneralDistress + longAnhedonia + longFears + site + sex + meds + race + inc + ethnicity')
        mdlco2 = fitlm(alldata,'CinguloOpercular  ~ longGeneralDistress + longAnhedonia + longFears  + site + sex + meds + race + inc + ethnicity')
    end

    %% AIM 2
    if mediation_networks == 1 
        % output doesn't save super nicely so I'm not going to bother. Will
        % run each line manually
        save_all_mdls = 0;
        % only writing code for VS connections to minimize tests
        if save_all_mdls == 1
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.Reward,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:RewardNetwork'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.Reward,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:RewardNetwork'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.Reward,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:RewardNetwork'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.DefaultMode,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:DefaultMode'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.DefaultMode,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:DefaultMode'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.DefaultMode,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:DefaultMode'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.FrontoParietal,'names',{'X:inflammation' 'Y:GeneralDistress' 'M:FrontoParietal'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.FrontoParietal,'names',{'X:inflammation' 'Y:Anhedonia' 'M:FrontoParietal'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.FrontoParietal,'names',{'X:inflammation' 'Y:Fears' 'M:FrontoParietal'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            
            resultsgendis = mediation(alldata.inflammation, alldata.longGeneralDistress, alldata.Salience,'names',{'X:inflammation' 'Y:GeneralDistress' 'M:Salience'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation, alldata.longAnhedonia, alldata.Salience,'names',{'X:inflammation' 'Y:Anhedonia' 'M:Salience'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation, alldata.longFears, alldata.Salience,'names',{'X:inflammation' 'Y:Fears' 'M:Salience'},'verbose','plots','doCIs','covs',[alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
        
        else
            keyboard
        end
    end

    % AIM 3 MODERATED MEDIATION
    if moderated_mediation_networks == 1
        % output doesn't save super nicely so I'm not going to bother. Will
        % run each line manually
        save_all_mdls = 0;
        % only writing code for VS connections to minimize tests
        if save_all_mdls == 1
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.Reward.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:Reward*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.Reward.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:Reward*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.Reward.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:Reward*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.DefaultMode.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:DefaultMode*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.DefaultMode.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:DefaultMode*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.DefaultMode.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:longitudinalsymptom' 'M:DefaultMode*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
                
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.FrontoParietal.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:GeneralDistress' 'M:FrontoParietal*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.FrontoParietal.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Anhedonia' 'M:FrontoParietal*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.FrontoParietal.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Fears' 'M:FrontoParietal*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
        
            resultsgendis = mediation(alldata.inflammation.*alldata.cti, alldata.longGeneralDistress, alldata.Salience.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:GeneralDistress' 'M:Salience*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsanhed = mediation(alldata.inflammation.*alldata.cti, alldata.longAnhedonia, alldata.Salience.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Anhedonia' 'M:Salience*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    
            resultsfears = mediation(alldata.inflammation.*alldata.cti, alldata.longFears, alldata.Salience.*alldata.cti,'names',{'X:inflammation*CTI' 'Y:Fears' 'M:Salience*CTI'},'verbose','plots','doCIs','covs',[alldata.cti,alldata.site,alldata.sex,alldata.ethnicity,alldata.race,alldata.inc,alldata.meds]);    

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
    check_networks = 0;
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
    listofseeds = [246 247];
    listoftargets = [244 245 ...
            105 108 111 116 117 118 ...
            102 110 122 204 206 208];
    
    if check_networks == 1
        atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
        % vs seitz idx = 246 247
        % amyg seitz idx = 244 245
        % ofc seitz idx = 105 111 116 117 118
        % acc seitz idx = 102 110 108 122 204 206 208
        listofregions = [246 247 244 245 ...
            105 111 116 117 118 ...
            102 108 110 122 204 206 208];
        hypatl = atl; hypatl.dat(:,1)=0;
        hypatlall = atl; hypatlall.dat(:,1)=0;
        for r = 1:length(listofregions)
            hypatl.dat(atl.dat==listofregions(r)) = 1;
            hypatlall.dat(atl.dat==listofregions(r)) = listofregions(r);
        end
    end
    %% important to change!!
    
    fnames(isnan(regressors.inflammation))=[];

    for f = 1:length(fnames)

        load(fnames{f});
        % this will load three variables. corr_vs_to_wholebrain, hypmat,
        % temp_brain. the one you want for this set of analyses is hypmat
        % which will be a VS voxel x 300 ROI masked connectivity matrix.
        % Available variables are: final_corr_acc_voxel_to_region
        % final_corr_ofc_voxel_to_region final_corr_amyg_voxel_to_region
        % final_corr_vs_voxel_to_region. For rest, corr_vs_voxel_to_region
        % and corresponding names for other seed regions. However, I'm
        % going in a different direction. Above, I confirmed the Seitzman
        % regions that correspond with all my regions of interest for this
        % project. I'm going to do a VS voxel x all other voxel corr matrix
        % in hyperalignment. This first round of hyperalignment will test
        % the idea that inflammation does have specific effects on
        % cortico-striatal-limbic connectivity, though it doesn't manifest
        % as an average (measured by standard metrics) but rather results
        % in subtle encoding differences at the voxel level. An alternative
        % I may want to try will be a voxelxtime hyperalignment.
        
        if rest == 1
            seeds = [seitz(246).all_data';seitz(247).all_data'];
            targets = [seitz(244).all_data';seitz(245).all_data';seitz(105).all_data';...
                seitz(111).all_data';seitz(116).all_data';seitz(117).all_data';...
                seitz(118).all_data';seitz(102).all_data';seitz(108).all_data';seitz(110).all_data';seitz(122).all_data';...
                seitz(204).all_data';seitz(206).all_data';seitz(208).all_data'];
            unaligned_mats{f} = atanh(corr(seeds',targets')');
        end
        clear temp
        if mid == 1
            seeds = [[seitz(1,246).all_data',seitz(1,246).all_data'];[seitz(1,247).all_data',seitz(1,247).all_data']];
            targets = [[seitz(1,244).all_data',seitz(2,244).all_data'];...
                [seitz(1,245).all_data',seitz(2,245).all_data'];...
                [seitz(1,105).all_data',seitz(2,105).all_data'];...
                [seitz(1,111).all_data',seitz(2,111).all_data'];...
                [seitz(1,116).all_data',seitz(2,116).all_data'];...
                [seitz(1,117).all_data',seitz(2,117).all_data'];...
                [seitz(1,118).all_data',seitz(2,118).all_data'];...
                [seitz(1,102).all_data',seitz(2,102).all_data'];...
                [seitz(1,108).all_data',seitz(2,108).all_data'];...
                [seitz(1,110).all_data',seitz(2,110).all_data'];...
                [seitz(1,122).all_data',seitz(2,122).all_data'];...
                [seitz(1,204).all_data',seitz(2,204).all_data'];...
                [seitz(1,206).all_data',seitz(2,206).all_data'];...
                [seitz(1,208).all_data',seitz(2,208).all_data']];
            unaligned_mats{f} = atanh(corr(seeds',targets')');
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
        curr_filename1 = fullfile(basedir,'/hyperalignment_results/',strcat('vs_data_perm_',num2str(perm),'.mat'));
        curr_filename2 = fullfile(basedir,'/hyperalignment_results/',strcat('vs_transforms_perm_',num2str(perm),'.mat'));
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
    if do_pls_regress == 1

        % easier to define an outcome here which will then be plugged in
        % throughout the document
        important_to_change = 1;
        outcome = regressors.longGeneralDistress(:,1);
        symptom = regressors.longGeneralDistress(:,1);
        symptom(isnan(regressors.inflammation))=[];
        outcome(isnan(regressors.inflammation))=[];
       
        for perm = 1:100
            load(analyze_fnames{perm})
            for sub = 1:length(aligned_data_resorted)
                final_aligned_data(:,sub) = aligned_data_resorted{sub}(:);
                final_unaligned_data(:,sub) = unaligned_mats{sub}(:); 
            end
            
            fprintf(strcat('Im on permutation: ',num2str(perm),'\n'))
            [aXL{perm},aYL{perm},aXS{perm},aYS{perm},aBETA{perm},aPCTVAR{perm},aMSE{perm},~] = plsregress(final_aligned_data',outcome,10,'CV',10);    
            clear final_aligned_data  
        end   

        for sub = 1:length(aligned_data_resorted)
            final_unaligned_data(:,sub) = unaligned_mats{sub}(:); 
        end

        [uXL,uYL,uXS,uYS,uBETA,uPCTVAR,uMSE,~] = plsregress(final_unaligned_data',outcome,10,'CV',10);
        % generate null distribution
        for i=1:100
            idx = randperm(size(outcome,1));
            fprintf(strcat('Im on permutation: ',num2str(i),'\n'))
            null_outcome = outcome(idx);
            [dXL{i},dYL{i},dXS{i},dYS{i},dBETA{i},dPCTVAR{i},dMSE{i},~] = plsregress(final_unaligned_data',null_outcome,10,'CV',10);                   
        end
        
        cd('shortcuts')
        save null_dist_model_results.mat dXL dYL dXS dYS dBETA dPCTVAR dMSE
        save aligned_dist_model_results.mat aXL aYL aXS aYS aBETA aPCTVAR aMSE
        save unaligned_dist_model_results.mat uXL uYL uXS uYS uBETA uPCTVAR uMSE
       
    end
end

if hyper_analyze_transforms == 1
    if mid == 1
        transform_files = filenames(fullfile('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/hyperalignment_results/*transforms*.mat'));
    end
    if rest == 1
        transform_files = filenames(fullfile('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices/hyperalignment_results/*transforms*.mat'));
    end

    important_to_change = 1;
    outcome = regressors.longGeneralDistress(:,1);
    symptom = regressors.longGeneralDistress(:,1);
    symptom(isnan(regressors.inflammation))=[];
    outcome(isnan(regressors.inflammation))=[];

    for perm = 1:length(transform_files)
        load(transform_files{perm})
        fprintf(strcat('Im doing perm: ', num2str(perm),'\n'))
        for sub = 1:length(transforms_resorted)

            % total magnitude combine rotation and translation
            totmagnitude_mat(:,:,sub) = transforms_resorted{sub}.T' * transforms_resorted{sub}.c';
            totmagnitude_clust(sub,:,:) = transforms_resorted{sub}.T' * transforms_resorted{sub}.c';
            totmagnitude_temp = (transforms_resorted{sub}.T' * transforms_resorted{sub}.c').^2;
            %totmagnitude_temp = (transforms_resorted{sub}.c * transforms_resorted{sub}.T).^2;
            totmagnitude_comp(sub,1,perm) = sum(totmagnitude_temp(:));
            totmagnitude_pls(sub,:,perm) = zscore(totmagnitude_temp(:));

            % magnitude of translation (distance metric)
            
            all_translation(sub,:,perm) = (sum((transforms_resorted{sub}.c.^2),2)).^0.5;     

            % voxelwise rotation (distance metric)
           
            all_rotation_mat(sub,:,perm) = zscore(mean(transforms_resorted{sub}.T.^2,2));
            % run model on voxelwise rotation
            

            % magnitude of scaling (single value)
            all_scale(sub,perm) = transforms_resorted{sub}.b;
            
            % determinant of rotation matrix (were data rotated with
            % reflection (det(T) = -1, or not?
            all_reflection(sub,perm) = det(transforms_resorted{sub}.T);

            
        end
        X = totmagnitude_pls(:,:,perm);
        
        [transformXL{perm},transformYL{perm},transformXS{perm},transformYS{perm},transformBETA{perm},transformPCTVAR{perm},transformMSE{perm},~] = plsregress(X,outcome,10,'CV',10);    
        
        idx = randperm(size(X,1));
        null_X = X(idx,:);
        [dtXL{perm},dtYL{perm},dtXS{perm},dtYS{perm},dtBETA{perm},dtPCTVAR{perm},dtMSE{perm},~] = plsregress(null_X,outcome,10,'CV',10);                   
        
        clear X  
        
        %tempfname = strcat('data_for_transform_mediation_and_cluster_perm',num2str(perm),'.mat');
        %save(fullfile(basedir,tempfname),"totmagnitude_clust")
    end
    
    % analyze transformation matrices. first total translation and rotation
    mdl_input = regressors;
    mdl_input(isnan(mdl_input.inflammation),:)=[];
    
    
    hyper_mdl_inputs = zscore([mean(all_translation(:,1),3),mean(all_scale(:,:),2),sum(all_reflection(:,:),2),mean(totmagnitude_comp,3)]);
    hyper_mdl_inputs = array2table(hyper_mdl_inputs); hyper_mdl_inputs.Properties.VariableNames = {'translation','scale','reflection','transrot'};
    mdl_input = [mdl_input,hyper_mdl_inputs];
    
    mdl_inf_trans = fitlm(mdl_input,'inflammation ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_inf_scale = fitlm(mdl_input,'inflammation ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_inf_transrot = fitlm(mdl_input,'inflammation ~ transrot + race + ethnicity + sex + meds + inc + site');

    mdl_gd_trans = fitlm(mdl_input,'GeneralDistress ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_gd_scale = fitlm(mdl_input,'GeneralDistress ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_gd_scale = fitlm(mdl_input,'GeneralDistress ~ transrot + race + ethnicity + sex + meds + inc + site');    

    mdl_anh_trans = fitlm(mdl_input,'Anhedonia ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_anh_scale = fitlm(mdl_input,'Anhedonia ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_anh_ref = fitlm(mdl_input,'Anhedonia ~ transrot + race + ethnicity + sex + meds + inc + site');

    mdl_fear_trans = fitlm(mdl_input,'Fears ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_fear_scale = fitlm(mdl_input,'Fears ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_fear_ref = fitlm(mdl_input,'Fears ~ transrot + race + ethnicity + sex + meds + inc + site');

    mdl_lgd_trans = fitlm(mdl_input,'longGeneralDistress ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lgd_scale = fitlm(mdl_input,'longGeneralDistress ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_lgd_ref = fitlm(mdl_input,'longGeneralDistress ~ transrot + race + ethnicity + sex + meds + inc + site');

    mdl_lanh_trans = fitlm(mdl_input,'longAnhedonia ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lanh_scale = fitlm(mdl_input,'longAnhedonia ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_lanh_ref = fitlm(mdl_input,'longAnhedonia ~ transrot + race + ethnicity + sex + meds + inc + site');

    mdl_lfear_trans = fitlm(mdl_input,'longFears ~ translation + race + ethnicity + sex + meds + inc + site');
    mdl_lfear_scale = fitlm(mdl_input,'longFears ~ scale + race + ethnicity + sex + meds + inc + site');
    mdl_lfear_ref = fitlm(mdl_input,'longFears ~ transrot + race + ethnicity + sex + meds + inc + site');
    
    
    

%     for perm = 1:100
%         
%         X = all_rotation_mat(:,:,perm);
%         
%         fprintf(strcat('Im on permutation: ',num2str(perm),'\n'))
%         [transformXL{perm},transformYL{perm},transformXS{perm},transformYS{perm},transformBETA{perm},transformPCTVAR{perm},transformMSE{perm},~] = plsregress(X,outcome,10,'CV',10);    
%         clear X  
%     end
%     
    if mid == 1
        cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/hyperalignment_results/shortcuts')
    end
    if rest == 1
        cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices/hyperalignment_results/shortcuts')
    end
%     
    
    save transforms_predict_model_results.mat transformXL transformYL transformXS transformYS transformBETA transformPCTVAR transformMSE totmagnitude_mat
    save null_transforms_predict_model_results.mat dtXL dtYL dtXS dtYS dtBETA dtPCTVAR dtMSE totmagnitude_mat

    
%     % connection wise test
%     count = 1;
%     significant_connections = zeros(size(totmagnitude_mat,1),size(totmagnitude_mat,2));
%     for i = 1:size(totmagnitude_mat,1)
%         %percent_complete = count/(size(totmagnitude_mat,1)*size(totmagnitude_mat,2));
% 
%         fprintf(strcat('Percentage done: ', num2str((count/58000)*100),'\n'))
% 
%         for j = 1:size(totmagnitude_mat,2)
%             
%                         
%             tempmdl = fitlm([mdl_input.inflammation],zscore(reshape(totmagnitude_mat(i,j,:),[size(totmagnitude_mat,3),1])));
%     
%             trans_pvalues(i,j) = tempmdl.Coefficients.pValue(2);
%             count = count + 1;
%     
%             
%         end
%     end

end

if visualize_pls_results == 1
    if mid == 1
        cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/hyperalignment_results/shortcuts')
    end
    if rest == 1
        cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/rest_corr_matrices/hyperalignment_results/shortcuts')
    end
    load aligned_dist_model_results.mat
    load unaligned_dist_model_results.mat
    load null_dist_model_results.mat
    load transforms_predict_model_results.mat
    load null_transforms_predict_model_results.mat

    for i = 1:100
        allMSE(i,1) = aMSE{i}(2,1);
        allMSE(i,2) = uMSE(2,1);
        allMSE(i,3) = dMSE{i}(2,1);
        allMSE(i,4) = transformMSE{i}(2,1);
        allMSE(i,5) = dtMSE{i}(2,1);

        % first component that I will try to analyze
        allPCT(i,1) = aPCTVAR{i}(2,1);
        allPCT(i,3) = dPCTVAR{i}(2,1);
        allPCT(i,2) = uPCTVAR(2,1);
        allPCT(i,4) = transformPCTVAR{i}(2,1);
        allPCT(i,5) = dtPCTVAR{i}(2,1);

        % all coverage
        totPCT(i,:,1) = aPCTVAR{i}(2,:);
        totPCT(i,:,3) = dPCTVAR{i}(2,:);
        totPCT(i,:,2) = uPCTVAR(2,:);
        totPCT(i,:,4) = transformPCTVAR{i}(2,:);
        totPCT(i,:,5) = dtPCTVAR{i}(2,:);

        loadings_to_visualize(:,:,i) = reshape(transformXL{i}(:,1),[828,71]);
        loadings_to_visualize_null(:,:,i) = reshape(dtXL{i}(:,1),[828,71]);
    end

    % visualize MSE
    % no differences in MSE. Distributions look super similar
    figure(); histogram(allMSE(:,1)); hold on; histogram(allMSE(:,3)); hold on; xline(allMSE(1,2)); hold on; histogram(allMSE(:,4)); hold on; histogram(allMSE(:,5));
    figure(); histogram(allMSE(:,1)); hold on; xline(allMSE(1,2)); hold on; histogram(allMSE(:,4));
    p_aligned_pls_mse = sum(allMSE(:,1)>allMSE(1,2))./100; % 
    p_transforms_pls_mse = sum(allMSE(:,4)>allMSE(1,2))./100;

    % visualize pct variance covered
    % unaligned data is no better than the null dist. Aligned data is WAY
    % better than both in terms of percent variance accounted for
    % pctall = figure(); histogram(allPCT(:,1)); hold on; histogram(allPCT(:,3)); hold on; xline(allPCT(1,2)); hold on; histogram(allPCT(:,4)); hold on; histogram(allPCT(:,5));
    pctall = figure(); histogram(allPCT(:,1)); hold on; xline(allPCT(1,2)); hold on; histogram(allPCT(:,4)); 
    %savefig(pctall,'pctall.fig'); savefig(pctall,'pct_hyperalignment.fig');
    p_aligned_pls_pct = sum(allPCT(:,1)<allPCT(1,2))./100; % p < 0.01
    p_transform_pls_pct = sum(allPCT(:,4)<allPCT(1,2))./100; % p < 0.01

    % visualize pls loadings of first component that accounts for 50-90% of
    % the variance in inflammation
    final_pls_loadings = mean(loadings_to_visualize,3);
    final_pls_loadings_null = mean(loadings_to_visualize_null,3);

    load('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/10001_final_ant.mat')
    idx_amyg = size(seitz(1,244).all_data,2) + size(seitz(1,245).all_data,2);
    idx_ofc = size(seitz(1,105).all_data,2) + size(seitz(1,108).all_data,2) + size(seitz(1,111).all_data,2) + size(seitz(1,116).all_data,2) + size(seitz(1,117).all_data,2) + size(seitz(1,118).all_data,2);
    idx_acc = size(seitz(1,102).all_data,2) + size(seitz(1,110).all_data,2) + size(seitz(1,122).all_data,2) + size(seitz(1,204).all_data,2) + size(seitz(1,206).all_data,2) + size(seitz(1,208).all_data,2);

    figure(); heatmap(sum(final_pls_loadings,2)')
    figure(); heatmap(zscore(sum(final_pls_loadings,2)'))

    % vs seitz idx = 246 247
    % amyg seitz idx = 244 245
    % ofc seitz idx = 105 111 116 117 118
    % acc seitz idx = 102 108 110 122 204 206 208
    

    idx_amyg = size(seitz(1,244).all_data,2) + size(seitz(1,245).all_data,2);
    idx_ofc = size(seitz(1,105).all_data,2) + size(seitz(1,111).all_data,2) + size(seitz(1,116).all_data,2) + size(seitz(1,117).all_data,2) + size(seitz(1,118).all_data,2);
    idx_acc = size(seitz(1,102).all_data,2) + size(seitz(1,108).all_data,2) + size(seitz(1,110).all_data,2) + size(seitz(1,122).all_data,2) + size(seitz(1,204).all_data,2) + size(seitz(1,206).all_data,2) + size(seitz(1,208).all_data,2);

    atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    % vs seitz idx = 246 247
    % amyg seitz idx = 244 245
    % ofc seitz idx = 105 111 116 117 118
    % acc seitz idx = 102 108 110 122 204 206 208
    listofregions = [244 245 ...
        105 111 116 117 118 ...
        102 108 110 122 204 206 208];
    namesofregions = {'amyg1' 'amyg2' 'ofc1' 'ofc2' 'ofc3' 'ofc4' 'ofc5' 'acc1' 'acc2' 'acc3' 'acc4' 'acc5' 'acc6' 'acc7'};
    hypatl_std = atl; hypatl_std.dat(:,1)=0;
    hypatl_avg = atl; hypatl_avg.dat(:,1)=0;
    hypatl_std_null = atl; hypatl_std_null.dat(:,1)=0;
    hypatl_avg_null = atl; hypatl_avg_null.dat(:,1)=0;

    
    sum_final_pls_loadings = sum(final_pls_loadings,2);
    sum_final_pls_loadings_null = sum(final_pls_loadings_null,2); 

    starting_idx = 1;

    % disimilarity of loadings_to_visualize
    for i = 1:100
        temp = mean(loadings_to_visualize(:,:,i),2);
        temp_null = mean(loadings_to_visualize_null(:,:,i),2);
        disim(:,i) = temp(:);
        disim_null(:,i) = temp_null(:);
    end
    
    disim_diff_perm = corr(disim) - corr(disim_null);
    disim_diff_vox = corr(disim') - corr(disim_null');

    for r = 1:length(listofregions)
        avg_pls(r) = mean(sum_final_pls_loadings(starting_idx:starting_idx+size(seitz(1,listofregions(r)-1).all_data,2)));
        std_pls(r) = std(sum_final_pls_loadings(starting_idx:starting_idx+size(seitz(1,listofregions(r)-1).all_data,2)));
        hypatl_avg.dat(atl.dat==listofregions(r)) = avg_pls(r);
        hypatl_std.dat(atl.dat==listofregions(r)) = std_pls(r);

        avg_pls_null(r) = mean(sum_final_pls_loadings_null(starting_idx:starting_idx+size(seitz(1,listofregions(r)-1).all_data,2)));
        std_pls_null(r) = std(sum_final_pls_loadings_null(starting_idx:starting_idx+size(seitz(1,listofregions(r)-1).all_data,2)));
        hypatl_avg_null.dat(atl.dat==listofregions(r)) = avg_pls_null(r);
        hypatl_std_null.dat(atl.dat==listofregions(r)) = std_pls_null(r);
        starting_idx = starting_idx + size(seitz(1,listofregions(r)).all_data,2);
    end
    
    %hypatl.mean         
    avg_z_table = avg_pls; avg_z_table=array2table(avg_z_table);
    avg_z_table.Properties.VariableNames = namesofregions;

    avg_z_table_null = avg_pls_null; avg_z_table_null=array2table(avg_z_table_null);
    avg_z_table_null.Properties.VariableNames = namesofregions;
     
    save avg_z_and_brain.mat avg_z_table hypatl_avg hypatl_std avg_z_table_null hypatl_avg_null hypatl_std_null avg_pls std_pls avg_pls_null std_pls_null loadings_to_visualize loadings_to_visualize_null
    
    
end

if mediation_with_hyp_transform == 1
    cd('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/mediation_and_clustering_data')
    load('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_corr_matrices/final_corr_mats/10001_final_ant.mat')
    med_fnames = filenames(fullfile('data*mat'));
    
    atl = fmri_data('/home/zaz3744/repo/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    
    important_to_change = 1;
    outcome = regressors.longGeneralDistress(:,1);
    symptom_long= regressors.longAnhedonia(:,1);
    symptom_long(isnan(regressors.inflammation))=[];
    symptom = regressors.Anhedonia(:,1);
    symptom(isnan(regressors.inflammation))=[];
    outcome(isnan(regressors.inflammation))=[];
    regressors(isnan(regressors.inflammation),:)=[];

    for perm = 1:length(med_fnames)
        load(med_fnames{perm});
        listofregions = [244 245 ...
        105 111 116 117 118 ...
        102 108 110 122 204 206 208];
        namesofregions = {'amyg1' 'amyg2' 'ofc1' 'ofc2' 'ofc3' 'ofc4' 'ofc5' 'acc1' 'acc2' 'acc3' 'acc4' 'acc5' 'acc6' 'acc7'};
        hypatl = atl; hypatl.dat(:,1)=0;
                
        starting_idx = 1;
        for r = 1:length(listofregions)
            
            current_data_temp{r} = mean(totmagnitude_clust(:,starting_idx:starting_idx+size(seitz(1,listofregions(r)).all_data,2)-1,:),3);
            starting_idx = starting_idx + size(seitz(1,listofregions(r)).all_data,2);
            
        end       
        totmagnitude_final(:,:,perm) = mean(totmagnitude_clust,3);
        bilateral_amyg(:,perm) = (mean(current_data_temp{1},2) + mean(current_data_temp{2},2))./2;
        bilateral_ofc(:,perm) = (mean(current_data_temp{3},2) +...
            mean(current_data_temp{4},2)+...
            mean(current_data_temp{5},2)+...
            mean(current_data_temp{6},2)+...
            mean(current_data_temp{7},2)...
            )./6;
        bilateral_acc(:,perm) = (mean(current_data_temp{8},2) +...
            mean(current_data_temp{9},2)+...
            mean(current_data_temp{10},2)+...
            mean(current_data_temp{11},2)+...
            mean(current_data_temp{12},2)+...
            mean(current_data_temp{13},2)+...
            mean(current_data_temp{14},2)...
            )./6;
       
    end
    results_amyg = mediation(outcome, symptom_long, zscore(mean(bilateral_amyg,2)),'M',symptom,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:Reward'},'verbose','plots','doCIs','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);    
    results_ofc = mediation(outcome, symptom_long, zscore(mean(bilateral_ofc,2)),'M',symptom,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:Reward'},'verbose','plots','doCIs','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);    
    results_acc = mediation(outcome, symptom_long, zscore(mean(bilateral_acc,2)),'M',symptom,'names',{'X:inflammation' 'Y:longitudinalsymptom' 'M:Reward'},'verbose','plots','doCIs','covs',[regressors.site,regressors.sex,regressors.ethnicity,regressors.race,regressors.inc,regressors.meds]);    

    [idx,C,sumd,D] = kmeans(zscore(mean(totmagnitude_final,3)),4);
    [idxclin,Cclin,sumdclin,Dclin] = kmeans([regressors.longGeneralDistress,regressors.longAnhedonia,regressors.longFears],4);

    idx1 = zeros(size(idx));idx1(idx==1)=1;
    idx2 = zeros(size(idx));idx2(idx==2)=1;
    idx3 = zeros(size(idx));idx3(idx==3)=1;
    idx4 = zeros(size(idx));idx4(idx==4)=1;


    [tbl,chi2,p,labels] = crosstab(idx,idxclin)

    clust_mdl = fitlm([idx,regressors.site,regressors.sex,regressors.meds,regressors.race,regressors.ethnicity,regressors.inc],regressors.longGeneralDistress)
    clust_mdl2 = fitlm(idx,regressors.longGeneralDistress)
    
end

