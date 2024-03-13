prep_confounds = 1;
prep_timing = 0;
%%
if prep_timing == 1
    savedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/final_betaseries_timing/';
    filedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/mid_spm_timing_baseline/';
    run = 'run-1';
    cd(fullfile(filedir,run))
    antfnames = filenames('anticipation/*.mat');
    confnames = filenames('consumption/*.mat');
    for files = 1:length(antfnames)
        anttemp = load(antfnames{files});
        contemp = load(confnames{files});
        pid = antfnames{files}(14:18);
        
        antevent = 'antgain';
        antindex = find(strcmp(anttemp.names(:),antevent));
        for i = 1:length(anttemp.onsets{antindex})
            durations1{1,i} = anttemp.durations{antindex}(i);
            names1{1,i} = anttemp.names{antindex};
            onsets1{1,i} = anttemp.onsets{antindex}(i);
        end

        conevent = 'hitcongain';
        conindex = find(strcmp(contemp.names(:),conevent));
        for j = 1:length(contemp.onsets{conindex})
            durations2{1,j} = contemp.durations{conindex}(j);
            names2{1,j} = contemp.names{conindex};
            onsets2{1,j} = contemp.onsets{conindex}(j);
        end
        
        durations = [durations1,durations2];
        onsets = [onsets1,onsets2];
        names = [names1,names2];
        

        filename = fullfile(savedir,strcat(pid,'_',run,'.mat'));
        save(filename,"durations","names","onsets");

        clear durations durations1 durations2 names names1 names2 onsets onsets1 onsets2 anttemp contemp
    end

end
%%
if prep_confounds == 1
    confounddir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/copy_fmriprep_confound';
    
    fnames = filenames('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/copy_fmriprep_confound/*rest*run-1*txt');
    
    ex1 =1; 
    ndummies = 2;
    for files = 1:length(fnames)
        T = readtable(fnames{files});
        pid = fnames{files}(78:82);
        outliers = table2array(T(:,contains(T.Properties.VariableNames,'motion')));
        transx = table2array(T(:,contains(T.Properties.VariableNames,'trans_x')));
        transy = table2array(T(:,contains(T.Properties.VariableNames,'trans_y')));
        transz = table2array(T(:,contains(T.Properties.VariableNames,'trans_z')));
        rotx = table2array(T(:,contains(T.Properties.VariableNames,'rot_x')));
        roty = table2array(T(:,contains(T.Properties.VariableNames,'rot_y')));
        rotz = table2array(T(:,contains(T.Properties.VariableNames,'rot_z')));
        gsr = table2array(T(:,contains(T.Properties.VariableNames,'global_signal')));
        csf = T.csf;
        wm = T.white_matter;
        fd(files) = nanmean(T.framewise_displacement);
    
        R = [transx,transy,transz,rotx,roty,rotz,gsr,outliers];
        R(isnan(R)) = 0;
        R = R(ndummies+1:size(R,1),:);
    
        if nanmean(T.framewise_displacement) > 0.25 
            pid_exclude_list{ex1,1} = pid;
            pid_exclude_list{ex1,2} = 'rest_run1';
            ex1 = ex1 + 1;
        end
        
        %save_name = fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/final_confound',strcat(pid,'_mid_run-1.mat'));
        %save(save_name,"R")
        
    end
end
