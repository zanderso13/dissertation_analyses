cd /Users/zacharyanderson/Documents/dissertation/phil_contrasts/ppi/

contrast = 'anticipation';
fnames1 = filenames(fullfile(strcat('sub*/ses-2/',contrast,'/run-1/con_0001.nii')));
fnames2 = filenames(fullfile(strcat('sub*/ses-2/',contrast,'/run-2/con_0001.nii')));

% find indices for run-1 and run-2 scans
for sub = 1:length(fnames2)
    PID{sub} = fnames2{sub}(1:9);
    ind(sub,:) = [sub,find(contains(fnames1,PID{sub}))];
end

exclusions = {'10001','10084','10094','10102','10125','10140','10143','10148','10161','10264','10274','10296','10319',...
  '10327','10422','10423','10434','10443','10461','10471',...
  '20032','20050','20108','20309','20317','20507','20564',...
  '20674','21111','21178','21223','21597','20877',... %,... % 20877 is missing specific trial types in one of the runs. 10272 only has one run
  '20384','21384',... % these two scans are duplicates so one of them has to be wrong. Excluding both for now.
  'sub-10006'	'sub-10074'	'sub-10111'	'sub-10141'	'sub-10161'	'sub-10196'	'sub-10236'	'sub-10264'	'sub-10272'	'sub-10274'	'sub-10282'	'sub-10296'	'sub-10308'	'sub-10341'	'sub-10434'	'sub-10459'	'sub-10481'	'sub-20150'	'sub-20309'	'sub-21597'	'sub-21684',... % motion exclusions >0.3
  'sub-10074'	'sub-10111'	'sub-10141'	'sub-10161'	'sub-10236'	'sub-10264'	'sub-10296'	'sub-10341'	'sub-20150'}; % motion exclusion in run 1 >0.3
  %'10272','10319',...
  %'10196','10326','21127','20556','10177','21384','20627','21386','10236','21238'}; % these exclusions are from the outlier detection

  
  
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
% plot(dfinal)

% now we'll take a look at framewise displacement





