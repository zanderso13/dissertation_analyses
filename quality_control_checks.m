% this is getting converted into a connectivity matrix script

cd /Users/zacharyanderson/Documents/ACNlab/BrainMAPD/PPI/dissertation_final_activation_data

fnames=;

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
  
excludethis = zeros(length(PID),1);
for sub = 1:length(exclusions)
    excludethistemp(:,1) = contains(PID,exclusions{sub});
    excludethis = excludethis + excludethistemp;
end
% apply exclusions
ind(excludethis>0,:) = [];
fnames2(excludethis>0,:) = [];
% get average mid information across runs
d1 = fmri_data(fnames1{ind(1,1)});
d2 = fmri_data(fnames1{ind(1,2)});
dfinal = d1; dfinal.dat = (dfinal.dat + d2.dat) ./ 2;
for sub = 2:length(ind)
    d1 = fmri_data(fnames1{ind(sub,1)});
    d2 = fmri_data(fnames1{ind(sub,2)});
    dfinal.dat(:,sub) = (d1.dat + d2.dat) ./ 2;
end

% this command checks the quality of the data
exclusions_from_plot = plot(dfinal);

% for plotting with montage
dfinal.X = ones(size(dfinal.dat,2),1);
reg_out = regress(dfinal); thresh_out = threshold(reg_out.t,0.05,'fdr','k',5);

% region of interest extraction. apply an atlas and check for outliers in
% connectivity/activation within a particular ROI

regions_to_extract = [157,158];
bilateral=1;
atl = fmri_data(fullfile('/Users/zacharyanderson/Documents/GitHub/SchizConnect/AAL3/AAL3v1.nii'));

if ~isempty(regions_to_extract)
    % create unilateral atlas from AAL
    atl_temp = atl;
    atl_temp.dat(atl_temp.dat>0) = 0;
    new_region_index = 1;
    for region = 1:length(regions_to_extract)
        atl_temp.dat(atl.dat==regions_to_extract(region)) = new_region_index;
        new_region_index = new_region_index + 1;
    end
    uni_atl = atl_temp; atl_temp = [];
    
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
    if bilateral == 1
        bi_atl.dat(bi_atl.dat > 0) = 1;
    end
end

r = extract_roi_averages(dfinal,bi_atl);

roi_quality_check = zscore(r.dat);

