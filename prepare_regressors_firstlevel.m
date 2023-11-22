prep_confounds = 0;
prep_timing = 1;
%%
if prep_timing == 1
    savedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/final_betaseries_timing/';
    filedir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/mid_spm_timing_baseline/';
    run = 'run-2';
    cd(fullfile(filedir,run))
    antfnames = filenames('anticipation/*.mat');
    confnames = filenames('consumption/*.mat');
    for files = 1:length(antfnames)
        anttemp = load(antfnames{files});
        contemp = load(confnames{files});
        pid = antfnames{files}(14:18);
        
        for events = 1:length(anttemp.names)-1
            for i = 1:length(anttemp.onsets{events})
                durations{i,events} = anttemp.durations{events}(i);
                names{i,events} = anttemp.names{events};
                onsets{i,events} = anttemp.onsets{events}(i);
                
            end

            for j = 1:length(contemp.onsets{events})
                durations2{j,events} = contemp.durations{events}(j);
                names2{j,events} = contemp.names{events};
                onsets2{j,events} = contemp.onsets{events}(j);
            end
        end
        
        durations = [durations(:)',durations2(:)'];
        names = [names(:)',names2(:)'];
        onsets = [onsets(:)',onsets2(:)'];

        durations = durations(~cellfun('isempty',durations));
        names = names(~cellfun('isempty',names));
        onsets = onsets(~cellfun('isempty',onsets));
        
        antdur = durations(strcmp(names, 'antgain')); 
        anton = onsets(strcmp(names, 'antgain')); 
        antname = names(strcmp(names, 'antgain'));
        condur = durations(strcmp(names, 'hitcongain'));
        conon = onsets(strcmp(names, 'hitcongain'));
        conname = names(strcmp(names, 'hitcongain'));
        
        clear durations names onsets

        durations = [antdur,condur];
        onsets = [anton,conon];
        names = [antname,conname];

        filename = fullfile(savedir,strcat(pid,'_',run,'.mat'));
        save(filename,"durations","names","onsets");

        clear durations names onsets antname conname anton conon antdur condur
    end

end
%%
if prep_confounds == 1
    confounddir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/copy_fmriprep_confound';
    
    fnames = filenames('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/copy_fmriprep_confound/*mid_run-2*txt');
    ex1 = 52; 
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
        csf = T.csf;
        wm = T.white_matter;
        fd(files) = nanmean(T.framewise_displacement);
    
        R = [transx,transy,transz,rotx,roty,rotz,csf,wm,outliers];
        R(isnan(R)) = 0;
        R = R(ndummies+1:size(R,1),:);
    
        if nanmean(T.framewise_displacement) > 0.2 
            pid_exclude_list{ex1,1} = pid;
            pid_exclude_list{ex1,2} = 'rest_run-1';
            ex1 = ex1 + 1;
        end
        save_name = fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/final_confound',strcat(pid,'_mid_run-2.mat'));
        save(save_name,"R")
        
    end
end
