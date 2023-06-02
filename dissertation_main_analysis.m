
load('/Users/zacharyanderson/Documents/dissertation/final_data/ofc_seed_mid_avg_data.mat');
load('/Users/zacharyanderson/Documents/dissertation/final_data/final_mid_fnames.mat');
load('/Users/zacharyanderson/Documents/dissertation/outcomes/demographics.mat');
load('/Users/zacharyanderson/Documents/dissertation/outcomes/trilevel.mat');
load('/Users/zacharyanderson/Documents/dissertation/outcomes/immune.mat');
load('/Users/zacharyanderson/Documents/dissertation/outcomes/meds.mat');

for sub = 1:length(fnames2)
    id{sub,1} = str2num(fnames2{sub}(5:9));
    if ~isempty(find(immune.PID==id{sub,1}))
        i(sub,1) = immune.T1BDicsavg(find(immune.PID==id{sub,1}));
    else
        fprintf(strcat(num2str(id{sub,1}),': missing immune data \n'))
        i(sub,1) = NaN;
    end
    if ~isempty(find(dem.PID==id{sub,1}))
        d(sub,1) = dem.sex(find(dem.PID==id{sub,1}));
    else
        fprintf(strcat(num2str(id{sub,1}),': missing dem data \n'))
        d(sub,1) = NaN;
    end
    if ~isempty(find(meds.PID==id{sub,1}))
        m(sub,1) = meds.T2SCcmidopany(find(meds.PID==id{sub,1}));
    else
        fprintf(strcat(num2str(id{sub,1}),': missing med data \n'))
        m(sub,1) = 0;
    end
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

site = cell2mat(id); site(site<20000) = 0; site(site>20000) = 1; site = array2table(site);
s = array2table(s); s.Properties.VariableNames = {'GeneralDistress','AnhedoniaApprehension','Fears'};
m = array2table(m); m.Properties.VariableNames = {'medication'};
i = array2table(i); i.Properties.VariableNames = {'inflammation'};
d = array2table(d); d.Properties.VariableNames = {'sex'};

final_outcomes = [s,m,i,d,site];

regressors_variablenames = {'Inflammation','GeneralDistress','AnhedoniaApprehension','Fears','sex','site','medication'};
regressors = [final_outcomes.inflammation,final_outcomes.GeneralDistress, final_outcomes.AnhedoniaApprehension, final_outcomes.Fears, final_outcomes.sex,final_outcomes.site,final_outcomes.medication];

% whole brain
dfinal.dat(:,isnan(regressors(:,1))|isnan(regressors(:,2))) = [];
regressors(isnan(regressors(:,1))|isnan(regressors(:,2)),:) = [];

dfinal.X = regressors;

statobj = regress(dfinal);
threshobj = threshold(statobj.t,0.001,'unc','k',10);

% orthviews(select_one_image(threshobj,2))

%% seed to seed

% load in AAL3 atlas and pull out what regions you want for targets
regions_to_extract = [157:158]; % pregenual ACC: 153:154, amygdala: 45:46, medial orb: 21:22
% caudate: 75:76, putamen: 77:78, NAcc: 157:158

atl = fmri_data(fullfile('/Users/zacharyanderson/Documents/GitHub/SchizConnect/AAL3/AAL3v1.nii'));
%atl = fmri_data(fullfile('/Users/zacharyanderson/Documents/dissertation/300_ROI_Set/ROIs_300inVol_MNI.nii'));
if ~isempty(regions_to_extract)  
    % create bilateral atlas from AAL
    atl_temp = atl;
    atl_temp.dat(atl_temp.dat>0) = 0;
    new_region_index = 1;
    for region = 1:2:length(regions_to_extract)
        atl_temp.dat(atl.dat==regions_to_extract(region)) = new_region_index;
        atl_temp.dat(atl.dat==regions_to_extract(region + 1)) = new_region_index;
        if sum(atl_temp.dat==new_region_index) == 0
            keyboard
        end
        new_region_index = new_region_index + 1;
    end
    bi_atl = atl_temp; atl_temp = [];
    
    bi_atl.dat(bi_atl.dat > 0) = 1;
    
end
vs = fmri_data('/Users/zacharyanderson/Documents/ACNlab/masks/VS_8mmsphere_Oldham_Rew.nii');
amyg = fmri_data('/Users/zacharyanderson/Documents/ACNlab/masks/HO_Amygdala_50prob.nii');
ofc = fmri_data('/Users/zacharyanderson/Documents/ACNlab/masks/OFC_8mmsphere_Oldham.nii');

r = extract_roi_averages(dfinal,vs);

dat_s2s = array2table(regressors); dat_s2s.Properties.VariableNames = regressors_variablenames;
target = r.dat; target = array2table(target); target.Properties.VariableNames = {'target_region'}; dat_s2s = [target,dat_s2s];


dat_s2s(isoutlier(dat_s2s.target_region,"gesd"),:) = [];
dat_s2s(isoutlier(dat_s2s.Inflammation,"gesd"),:) = [];
dat_s2s(isoutlier(dat_s2s.GeneralDistress,"gesd"),:) = [];
dat_s2s(isoutlier(dat_s2s.AnhedoniaApprehension,"gesd"),:) = [];
dat_s2s(isoutlier(dat_s2s.Fears,"gesd"),:) = [];

% whole sample
fitlm(dat_s2s,'target_region ~ Inflammation + medication + site')
fitlm(dat_s2s,'target_region ~ Inflammation*GeneralDistress + medication + site')
fitlm(dat_s2s,'target_region ~ GeneralDistress*Inflammation + AnhedoniaApprehension + Fears + medication + site') 


% let's separate male/female assigned at birth

dat_s2s1 = dat_s2s; dat_s2s1(dat_s2s1.sex==1,:) = []; % assigned male at birth
dat_s2s2 = dat_s2s; dat_s2s2(dat_s2s2.sex==0,:) = []; % assigned female at birth

fitlm(dat_s2s1,'target_region ~ GeneralDistress + AnhedoniaApprehension + Fears + medication + site') 
fitlm(dat_s2s2,'target_region ~ Inflammation + medication + site') 

% work shows it's a subset of individuals with mood disorders who are
% impacted by inflammation related effects on neural circuits. I will try
% to break that down below with a series of interaction models

fitlm(dat_s2s1,'target_region ~ Inflammation*GeneralDistress + medication + site') % male at birth
fitlm(dat_s2s2,'target_region ~ Inflammation*GeneralDistress + medication + site') % female at birth

fitlm(dat_s2s1,'target_region ~ Inflammation*AnhedoniaApprehension + medication + site') % male at birth
fitlm(dat_s2s2,'target_region ~ Inflammation*AnhedoniaApprehension + medication + site') % female at birth

fitlm(dat_s2s1,'target_region ~ Inflammation*Fears + medication + site') % male at birth
fitlm(dat_s2s2,'target_region ~ Inflammation*Fears + medication + site') % female at birth



