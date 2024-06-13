
mid = 1;
rest = 0;
visualize_pls_results = 0;
compare_loadings_across_outcomes = 0;
plot_voxelwise_transforms_with_histograms_of_corr = 0;
rsa_voxelwise_transform_correlations = 1;
do_heatmaps_for_whole_networks = 0;

control_visual = 1;

simulations = 1;

display_regressors = 0;

if visualize_pls_results == 1
    if mid == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/anhedonia')
    end
    if rest == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/rest_shortcuts/anhedonia')
    end
    load aligned_dist_model_results.mat
    load unaligned_dist_model_results.mat
    load null_dist_model_results.mat
    load transforms_predict_model_results.mat

    for i = 1:100
        allMSE(i,1) = aMSE{i}(2,1);
        allMSE(i,2) = uMSE(2,1);
        allMSE(i,3) = dMSE{i}(2,1);
        allMSE(i,4) = transformMSE{i}(2,1);

        % first component that I will try to analyze
        allPCT(i,1) = aPCTVAR{i}(2,1);
        allPCT(i,3) = dPCTVAR{i}(2,1);
        allPCT(i,2) = uPCTVAR(2,1);
        allPCT(i,4) = transformPCTVAR{i}(2,1);

        % all coverage
        totPCT(i,:,1) = aPCTVAR{i}(2,:);
        totPCT(i,:,3) = dPCTVAR{i}(2,:);
        totPCT(i,:,2) = uPCTVAR(2,:);
        totPCT(i,:,4) = transformPCTVAR{i}(2,:);
    end

    % visualize MSE
    % no differences in MSE. Distributions look super similar
    figure(); histogram(allMSE(:,1)); hold on; xline(allMSE(1,2),'LineWidth', 5); hold on; histogram(allMSE(:,4));
    figure(); histogram(allMSE(:,1)); hold on; xline(allMSE(1,2),'LineWidth', 5); hold on; histogram(allMSE(:,4));
    p_aligned_pls_mse = sum(allMSE(:,1)>allMSE(1,2))./100; % 
    p_transforms_pls_mse = sum(allMSE(:,4)>allMSE(1,2))./100;

    % visualize pct variance covered
    % unaligned data is no better than the null dist. Aligned data is WAY
    % better than both in terms of percent variance accounted for
    pctall = figure(); histogram(allPCT(:,1)); hold on; xline(allPCT(1,2),'LineWidth', 5); hold on; histogram(allPCT(:,4)); 
    pcthyp = figure(); histogram(allPCT(:,1)); hold on; xline(allPCT(1,2),'LineWidth', 5); hold on; histogram(allPCT(:,4));
    savefig(pctall,'pctall.fig'); savefig(pcthyp,'pct_hyperalignment.fig');
    p_aligned_pls_pct = sum(allPCT(:,1)<allPCT(1,2))./100; % p < 0.01
    p_transform_pls_pct = sum(allPCT(:,4)<allPCT(1,2))./100; % p < 0.01
    
    
    load('avg_z_and_brain.mat')
    %figure(); montage(hypatl_avg)
    %figure(); montage(hypatl_std)

    if mid == 1
        hypatl_avg.fullpath = '/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/anhedonia/hypatl_avg.nii';
        write(hypatl_avg,'overwrite')
    end
    if rest == 1
        hypatl_avg.fullpath = '/Users/zacharyanderson/Documents/ACNlab/dissertation_final/rest_shortcuts/anhedonia/hypatl_avg.nii';
        write(hypatl_avg,'overwrite')    
    end
       

end

if compare_loadings_across_outcomes == 1
    fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/*shortcuts/*/avg_z_and_brain.mat'));
    amid = load(fnames{1});
    gmid = load(fnames{2});
    imid = load(fnames{3});
    arest = load(fnames{4});
    grest = load(fnames{5});
    irest = load(fnames{6});
    
    sim_mat = [amid.avg_z_table;gmid.avg_z_table;imid.avg_z_table;arest.avg_z_table;grest.avg_z_table;irest.avg_z_table];
    sim_mat.Properties.RowNames = {'Anhedonia MID','GD MID','Inf MID','Anhedonia Rest','GD rest','Inf rest'};
    
    % similarity of features
    figure();
    featfig = heatmap(corr(table2array(sim_mat))); featfig.XDisplayLabels = sim_mat.Properties.VariableNames;
    featfig.YDisplayLabels = sim_mat.Properties.VariableNames;
    % similarity of outcomes
    figure();
    outfig = heatmap(corr(table2array(sim_mat)')); outfig.XDisplayLabels = sim_mat.Properties.RowNames;
    outfig.YDisplayLabels = sim_mat.Properties.RowNames;

end


if plot_voxelwise_transforms_with_histograms_of_corr == 1

    if mid == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/anhedonia')
    end
    if rest == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/rest_shortcuts/anhedonia')
    end
    
    load('corr_matrices_for_heatmaps.mat');
    load('transforms_to_plot.mat');

    atl = fmri_data('/Users/zacharyanderson/Documents/GitHub/dissertation_analyses/300_ROI_Set/ROIs_300inVol_MNI.nii');
    listofregions = [244 245 ...
        105 111 116 117 118 ...
        102 108 110 122 204 206 208];
    namesofregions = {'amyg1' 'amyg2' 'ofc1' 'ofc2' 'ofc3' 'ofc4' 'ofc5' 'acc1' 'acc2' 'acc3' 'acc4' 'acc5' 'acc6' 'acc7'};
    
    figure();
    subplot(1,2,1);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'amyg1'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='Amygdala 1';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,2,2);
    h2 = heatmap(round(corr_transform{find(strcmp(namesofregions,'amyg2'))},3),'CellLabelColor','none','Colormap',parula);
    h2.XDisplayLabels='Amygdala 2';
    h2.YDisplayLabels = nan(size(h2.YDisplayData));

    figure();
    subplot(1,2,1);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'amyg1'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'amyg1'))}),3),'LineWidth',3)
    xlabel('Amygdala 1')
    subplot(1,2,2);
    hi2 = histogram(round(corr_transform{find(strcmp(namesofregions,'amyg2'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'amyg2'))}),3),'LineWidth',3)
    xlabel('Amygdala 2')

    atlamyg = atl; atlamyg.dat(:,1)=0; atlamyg.dat(atl.dat==244) = 1; atlamyg.dat(atl.dat==245) = 2;
    atlamyg.fullpath='/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/amyg.nii';
    write(atlamyg,'overwrite')

    figure();
    subplot(1,5,1);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'ofc1'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='mOFC 1';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,5,2);
    h2 = heatmap(round(corr_transform{find(strcmp(namesofregions,'ofc2'))},3),'CellLabelColor','none','Colormap',parula);
    h2.XDisplayLabels='mOFC 2';
    h2.YDisplayLabels = nan(size(h2.YDisplayData));
    subplot(1,5,3);
    h3 = heatmap(round(corr_transform{find(strcmp(namesofregions,'ofc3'))},3),'CellLabelColor','none','Colormap',parula);
    h3.XDisplayLabels='mOFC 3';
    h3.YDisplayLabels = nan(size(h3.YDisplayData));
    subplot(1,5,4);
    h4 = heatmap(round(corr_transform{find(strcmp(namesofregions,'ofc4'))},3),'CellLabelColor','none','Colormap',parula);
    h4.XDisplayLabels='mOFC 4';
    h4.YDisplayLabels = nan(size(h4.YDisplayData));

    %figure();
    subplot(1,5,5);
    h5 = heatmap(round(corr_transform{find(strcmp(namesofregions,'ofc5'))},3),'CellLabelColor','none','Colormap',parula);
    h5.XDisplayLabels='mOFC 5';
    h5.YDisplayLabels = nan(size(h5.YDisplayData));

    figure();
    subplot(1,5,1);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'ofc1'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'ofc1'))}),3),'LineWidth',3)
    xlabel('mOFC 1')
    subplot(1,5,2);
    hi2 = histogram(round(corr_transform{find(strcmp(namesofregions,'ofc2'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'ofc2'))}),3),'LineWidth',3)
    xlabel('mOFC 2')
    subplot(1,5,3);
    hi3 = histogram(round(corr_transform{find(strcmp(namesofregions,'ofc3'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'ofc3'))}),3),'LineWidth',3)
    xlabel('mOFC 3')
    subplot(1,5,4);
    hi4 = histogram(round(corr_transform{find(strcmp(namesofregions,'ofc4'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'ofc4'))}),3),'LineWidth',3)
    xlabel('mOFC 4')

    %figure();
    subplot(1,5,5);
    hi5 = histogram(round(corr_transform{find(strcmp(namesofregions,'ofc5'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'ofc5'))}),3),'LineWidth',3)
    xlabel('mOFC 5')

    atlofc = atl; atlofc.dat(:,1)=0; atlofc.dat(atl.dat==105) = 1; atlofc.dat(atl.dat==111) = 2;
    atlofc.dat(atl.dat==116) = 3; atlofc.dat(atl.dat==117) = 4; atlofc.dat(atl.dat==118) = 5;
    atlofc.fullpath='/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/ofc.nii';
    write(atlofc,'overwrite')

    figure();
    subplot(1,7,1);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc1'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 1';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,7,2)
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc2'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 2';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,7,3);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc3'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 3';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,7,4);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc4'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 4';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));

    %figure();
    subplot(1,7,5);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc5'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 5';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,7,6);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc6'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 6';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));
    subplot(1,7,7);
    h1 = heatmap(round(corr_transform{find(strcmp(namesofregions,'acc7'))},3),'CellLabelColor','none','Colormap',parula);
    h1.XDisplayLabels='ACC 7';
    h1.YDisplayLabels = nan(size(h1.YDisplayData));

    figure();
    subplot(1,7,1);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc1'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc1'))}),3),'LineWidth',3)
    xlabel('ACC 1')
    subplot(1,7,2);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc2'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc2'))}),3),'LineWidth',3)
    xlabel('ACC 2')
    subplot(1,7,3);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc3'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc3'))}),3),'LineWidth',3)
    xlabel('ACC 3')
    subplot(1,7,4);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc4'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc4'))}),3),'LineWidth',3)
    xlabel('ACC 4')

    %figure();
    subplot(1,7,5);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc5'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc5'))}),3),'LineWidth',3)
    xlabel('ACC 5')
    subplot(1,7,6);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc6'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc6'))}),3),'LineWidth',3)
    xlabel('ACC 6')
    subplot(1,7,7);
    hi1 = histogram(round(corr_transform{find(strcmp(namesofregions,'acc7'))},3)); hold on;
    xline(round(mean(corr_transform{find(strcmp(namesofregions,'acc7'))}),3),'LineWidth',3)
    xlabel('ACC 7')

    atlacc = atl; atlacc.dat(:,1)=0; atlacc.dat(atl.dat==102) = 1; atlacc.dat(atl.dat==108) = 2;
    atlacc.dat(atl.dat==110) = 3; atlacc.dat(atl.dat==122) = 4; atlacc.dat(atl.dat==204) = 5;
    atlacc.dat(atl.dat==206) = 6; atlacc.dat(atl.dat==208) = 7;
    atlacc.fullpath='/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/acc.nii';
    write(atlacc,'overwrite')
    
    table_to_display = [mean(corr_transform{find(strcmp(namesofregions,'amyg1'))}),std(corr_transform{find(strcmp(namesofregions,'amyg1'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'amyg2'))}),std(corr_transform{find(strcmp(namesofregions,'amyg2'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'ofc1'))}),std(corr_transform{find(strcmp(namesofregions,'ofc1'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'ofc2'))}),std(corr_transform{find(strcmp(namesofregions,'ofc2'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'ofc3'))}),std(corr_transform{find(strcmp(namesofregions,'ofc3'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'ofc4'))}),std(corr_transform{find(strcmp(namesofregions,'ofc4'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'ofc5'))}),std(corr_transform{find(strcmp(namesofregions,'ofc5'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc1'))}),std(corr_transform{find(strcmp(namesofregions,'acc1'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc2'))}),std(corr_transform{find(strcmp(namesofregions,'acc2'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc3'))}),std(corr_transform{find(strcmp(namesofregions,'acc3'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc4'))}),std(corr_transform{find(strcmp(namesofregions,'acc4'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc5'))}),std(corr_transform{find(strcmp(namesofregions,'acc5'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc6'))}),std(corr_transform{find(strcmp(namesofregions,'acc6'))}),...
    mean(corr_transform{find(strcmp(namesofregions,'acc7'))}),std(corr_transform{find(strcmp(namesofregions,'acc7'))})];

    table_to_display = array2table(table_to_display);
    table_to_display.Properties.VariableNames = {'mean amyg 1','sd amyg 1','mean amyg 2','sd amyg 2','mean ofc 1','sd ofc 1',...
        'mean ofc 2','sd ofc 2','mean ofc 3','sd ofc 3','mean ofc 4','sd ofc 4','mean ofc 5','sd ofc 5',...
        'mean acc 1','sd acc 1','mean acc 2','sd acc 2','mean acc 3','sd acc 3','mean acc 4','sd acc 4','mean acc 5','sd acc 5',...
        'mean acc 6','sd acc 6','mean acc 7','sd acc 7'};
end

if rsa_voxelwise_transform_correlations == 1
    
    cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final')    

    fnames = filenames(fullfile('*shortcuts*/*/corr_transform.mat'));
    namesofregions = {'amyg1' 'amyg2' 'ofc1' 'ofc2' 'ofc3' 'ofc4' 'ofc5' ...
        'acc1' 'acc2' 'acc3' 'acc4' 'acc5' 'acc6' 'acc7'};

    amid = load(fnames{1}); 
    gmid = load(fnames{2});
    imid = load(fnames{3});
    arest = load(fnames{4});
    grest = load(fnames{5});
    irest = load(fnames{6});

    for i=1:14
        amid.corr_transform{i} = amid.corr_transform{i}(:).*-1;
        arest.corr_transform{i} = arest.corr_transform{i}(:).*-1;
    end

    amyg1 = [amid.corr_transform{1}(:),gmid.corr_transform{1}(:),imid.corr_transform{1}(:),arest.corr_transform{1}(:),grest.corr_transform{1}(:),irest.corr_transform{1}(:)];
    amyg2 = [amid.corr_transform{2}(:),gmid.corr_transform{2}(:),imid.corr_transform{2}(:),arest.corr_transform{2}(:),grest.corr_transform{2}(:),irest.corr_transform{2}(:)];
    figure();
    amygall = [amyg1;amyg2];
    
    [bootstat_amyg_gi_mid,bootsam_amyg_gi_mid] = bootstrp(1000,@corr,amygall(:,2)',amygall(:,3)');
    [bootstat_amyg_ai_mid,bootsam_amyg_ai_mid] = bootstrp(1000,@corr,amygall(:,1)',amygall(:,3)');
    [bootstat_amyg_ag_mid,bootsam_amyg_ag_mid] = bootstrp(1000,@corr,amygall(:,1)',amygall(:,2)');
    [bootstat_amyg_gi_rest,bootsam_amyg_gi_rest] = bootstrp(1000,@corr,amygall(:,5)',amygall(:,6)');
    [bootstat_amyg_ai_rest,bootsam_amyg_ai_rest] = bootstrp(1000,@corr,amygall(:,4)',amygall(:,6)');
    [bootstat_amyg_ag_rest,bootsam_amyg_ag_rest] = bootstrp(1000,@corr,amygall(:,4)',amygall(:,5)');

    p_boot_amyg_ai_rest = sum(bootstat_amyg_ai_rest<0)/1000;
    p_boot_amyg_ai_mid = sum(bootstat_amyg_ai_mid>0)/1000;

    p_boot_amyg_gi_rest = sum(bootstat_amyg_gi_rest>0)/1000;
    p_boot_amyg_gi_mid = sum(bootstat_amyg_gi_mid>0)/1000;

    p_boot_amyg_ag_rest = sum(bootstat_amyg_ag_rest>0)/1000;
    p_boot_amyg_ag_mid = sum(bootstat_amyg_ag_mid>0)/1000;
    
    figure();
    amygh = heatmap(corr(amygall),'ColorMap',parula); amygh.XDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    amygh.YDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    amygh.Title = 'VS-Amygdala circuit transformations';
 
    ofc1 = [amid.corr_transform{3}(:),gmid.corr_transform{3}(:),imid.corr_transform{3}(:),arest.corr_transform{3}(:),grest.corr_transform{3}(:),irest.corr_transform{3}(:)];
    ofc2 = [amid.corr_transform{4}(:),gmid.corr_transform{4}(:),imid.corr_transform{4}(:),arest.corr_transform{4}(:),grest.corr_transform{4}(:),irest.corr_transform{4}(:)];
    ofc3 = [amid.corr_transform{5}(:),gmid.corr_transform{5}(:),imid.corr_transform{5}(:),arest.corr_transform{5}(:),grest.corr_transform{5}(:),irest.corr_transform{5}(:)];
    ofc4 = [amid.corr_transform{6}(:),gmid.corr_transform{6}(:),imid.corr_transform{6}(:),arest.corr_transform{6}(:),grest.corr_transform{6}(:),irest.corr_transform{6}(:)];
    ofc5 = [amid.corr_transform{7}(:),gmid.corr_transform{7}(:),imid.corr_transform{7}(:),arest.corr_transform{7}(:),grest.corr_transform{7}(:),irest.corr_transform{7}(:)];

    ofcall = [ofc1;ofc2;ofc3;ofc4;ofc5];

    [bootstat_ofc_gi_mid,bootsam_ofc_gi_mid] = bootstrp(1000,@corr,ofcall(:,2)',ofcall(:,3)');
    [bootstat_ofc_ai_mid,bootsam_ofc_ai_mid] = bootstrp(1000,@corr,ofcall(:,1)',ofcall(:,3)');
    [bootstat_ofc_gi_rest,bootsam_ofc_gi_rest] = bootstrp(1000,@corr,ofcall(:,5)',ofcall(:,6)');
    [bootstat_ofc_ai_rest,bootsam_ofc_ai_rest] = bootstrp(1000,@corr,ofcall(:,4)',ofcall(:,6)');
    

    p_boot_ofc_ai_rest = sum(bootstat_ofc_ai_rest<0)/1000;
    p_boot_ofc_ai_mid = sum(bootstat_ofc_ai_mid>0)/1000;

    p_boot_ofc_gi_rest = sum(bootstat_ofc_gi_rest>0)/1000;
    p_boot_ofc_gi_mid = sum(bootstat_ofc_gi_mid>0)/1000;

    figure();
    ofch = heatmap(corr(ofcall),'ColorMap',parula); ofch.XDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    ofch.YDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    ofch.Title = 'VS-mOFC circuit transformations';
   
    acc1 = [amid.corr_transform{8}(:),gmid.corr_transform{8}(:),imid.corr_transform{8}(:),arest.corr_transform{8}(:),grest.corr_transform{8}(:),irest.corr_transform{8}(:)];
    acc2 = [amid.corr_transform{9}(:),gmid.corr_transform{9}(:),imid.corr_transform{9}(:),arest.corr_transform{9}(:),grest.corr_transform{9}(:),irest.corr_transform{9}(:)];
    acc3 = [amid.corr_transform{10}(:),gmid.corr_transform{10}(:),imid.corr_transform{10}(:),arest.corr_transform{10}(:),grest.corr_transform{10}(:),irest.corr_transform{10}(:)];
    acc4 = [amid.corr_transform{11}(:),gmid.corr_transform{11}(:),imid.corr_transform{11}(:),arest.corr_transform{11}(:),grest.corr_transform{11}(:),irest.corr_transform{11}(:)];
    acc5 = [amid.corr_transform{12}(:),gmid.corr_transform{12}(:),imid.corr_transform{12}(:),arest.corr_transform{12}(:),grest.corr_transform{12}(:),irest.corr_transform{12}(:)];
    acc6 = [amid.corr_transform{13}(:),gmid.corr_transform{13}(:),imid.corr_transform{13}(:),arest.corr_transform{13}(:),grest.corr_transform{13}(:),irest.corr_transform{13}(:)];
    acc7 = [amid.corr_transform{14}(:),gmid.corr_transform{14}(:),imid.corr_transform{14}(:),arest.corr_transform{14}(:),grest.corr_transform{14}(:),irest.corr_transform{14}(:)];
    
    accall = [acc1;acc2;acc3;acc4;acc5;acc6;acc7];
    figure();
    acch = heatmap(corr(accall),'ColorMap',parula); acch.XDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    acch.YDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    acch.Title = 'VS-ACC circuit transformations';

    [bootstat_acc_gi_mid,bootsam_acc_gi_mid] = bootstrp(1000,@corr,accall(:,2)',accall(:,3)');
    [bootstat_acc_ai_mid,bootsam_acc_ai_mid] = bootstrp(1000,@corr,accall(:,1)',accall(:,3)');
    [bootstat_acc_gi_rest,bootsam_acc_gi_rest] = bootstrp(1000,@corr,accall(:,5)',accall(:,6)');
    [bootstat_acc_ai_rest,bootsam_acc_ai_rest] = bootstrp(1000,@corr,accall(:,4)',accall(:,6)');

    p_boot_acc_ai_rest = sum(bootstat_acc_ai_rest>0)/1000;
    p_boot_acc_ai_mid = sum(bootstat_acc_ai_mid>0)/1000;

    p_boot_acc_gi_rest = sum(bootstat_acc_gi_rest>0)/1000;
    p_boot_acc_gi_mid = sum(bootstat_acc_gi_mid>0)/1000;
    
    wholebrain = [amygall;ofcall;accall];

    figure();
    wholebrainh = heatmap(corr(wholebrain),'ColorMap',parula); wholebrainh.XDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};
    wholebrainh.Title = 'All cortico-striatal-amygdala circuit transformations';
    wholebrainh.YDisplayLabels={'Anhedonia MID','GD MID','Inflammation MID','Anhedonia Rest','GD Rest','Inflammation Rest'};

    figure();
    subplot(3,1,1)
    histogram(bootstat_amyg_gi_rest); hold on; histogram(bootstat_amyg_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-amygdala')
    %title('Models relating General Distress and Inflammation with rsfMRI data')
    subplot(3,1,2)
    histogram(bootstat_ofc_gi_rest); hold on; histogram(bootstat_ofc_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-mOFC')
    subplot(3,1,3)
    histogram(bootstat_acc_gi_rest); hold on; histogram(bootstat_acc_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-ACC')

    figure();
    subplot(3,1,1)
    histogram(bootstat_amyg_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-amygdala')
    %title('Models relating Anhedonia-Apprehension and Inflammation with rsfMRI data')
    subplot(3,1,2)
    histogram(bootstat_ofc_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-mOFC')
    subplot(3,1,3)
    histogram(bootstat_acc_ai_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-ACC')

    figure();
    subplot(3,1,1)
    histogram(bootstat_amyg_gi_mid); hold on; histogram(bootstat_amyg_ai_mid);
    xlabel('Correlation coefficients')
    ylabel('VS-amygdala')
    %title('Models relating Anhedonia-Apprehension and Inflammation with MID data')
    subplot(3,1,2)
    histogram(bootstat_ofc_gi_mid); hold on; histogram(bootstat_ofc_ai_mid);
    xlabel('Correlation coefficients')
    ylabel('VS-mOFC')
    subplot(3,1,3)
    histogram(bootstat_acc_gi_mid); hold on; histogram(bootstat_acc_ai_mid);
    xlabel('Correlation coefficients')
    ylabel('VS-ACC')
    
    % final histogram to assess for sig differences between GD and Anh
    % using VS-amygdala
    figure();
    subplot(2,1,1)
    histogram(bootstat_amyg_ag_rest);
    xlabel('Correlation coefficients')
    ylabel('VS-Amyg during resting state')
    subplot(2,1,2)
    histogram(bootstat_amyg_ag_mid);
    xlabel('Correlation coefficients')
    ylabel('VS-Amyg during MID task')

end

if do_heatmaps_for_whole_networks == 1
    mid_heat = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/whole_seitz_mat.mat');
    rest_heat = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/rest_shortcuts/whole_seitz_mat.mat');
    figure();
    heatmap(nanmean(mid_heat.final_mat,3),'ColorMap',parula)
    figure();
    heatmap(nanmean(rest_heat.final_mat,3),'ColorMap',parula)
    figure();
    heatmap(nanmean(mid_heat.final_mat,3)-nanmean(rest_heat.final_mat,3),'ColorMap',parula)
end

if simulations == 1

    load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/hyperalignment_simulation/pls_u.mat')
    load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/hyperalignment_simulation/pls_a.mat')
    load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/hyperalignment_simulation/pls_t.mat')
    load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/hyperalignment_simulation/rsa.mat')

    % 1 - inflammation
    % 2 - anhedonia
    % 3 - general distress

    for perm = 1:100
        pcta_inflam(perm,1) = aPCTVAR1{perm}(2,1);
        pcta_anhed(perm,1) = aPCTVAR2{perm}(2,1);
        pcta_gendis(perm,1) = aPCTVAR3{perm}(2,1);

        pctu_inflam(perm,1) = uPCTVAR1{perm}(2,1);
        pctu_anhed(perm,1) = uPCTVAR2{perm}(2,1);
        pctu_gendis(perm,1) = uPCTVAR3{perm}(2,1);

        pctt_inflam(perm,1) = tPCTVAR1{perm}(2,1);
        pctt_anhed(perm,1) = tPCTVAR2{perm}(2,1);
        pctt_gendis(perm,1) = tPCTVAR3{perm}(2,1);

        [bootstat_ai_sim(:,perm),~] = bootstrp(10,@corr,corr_map(:,1,perm),corr_map(:,2,perm));
        [bootstat_gi_sim(:,perm),~] = bootstrp(10,@corr,corr_map(:,1,perm),corr_map(:,3,perm));
        [bootstat_ag_sim(:,perm),~] = bootstrp(10,@corr,corr_map(:,2,perm),corr_map(:,3,perm));
    end
    
    figure()
    subplot(3,1,1);
    histogram(pctt_inflam);
    title('Transformations of random data predicting inflammation')
    subplot(3,1,2);
    histogram(pctt_gendis);
    title('Transformations of random data predicting General Distress')
    subplot(3,1,3);
    histogram(pctt_anhed);
    title('Transformations of random data predicting Anhedonia-Apprehension')


    figure()
    subplot(2,1,1);
    title('RSA: Anhedonia and Inflammation')
    histogram(bootstat_ai_sim(:));
    subplot(2,1,2);
    title('RSA: General Distress and Inflammation')
    histogram(bootstat_gi_sim(:));

    if rsa_voxelwise_transform_correlations==1
        p = 0.001;
        %%% Resting state
        %% ACC
        [p_gi_acc1,h_gi_acc1,stat_gi_acc1] = ranksum(bootstat_gi_sim(:), bootstat_acc_gi_rest,'Alpha',p);

        % Display results
        if h_gi_acc1 == 0
            fprintf('ACC comparing general distress and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_gi_acc1);
        else
            fprintf('ACC comparing general distress and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_gi_acc1);
        end

        [p_ai_acc1,h_ai_acc1,stat_ai_acc1] = ranksum(bootstat_ai_sim(:), bootstat_acc_ai_rest,'Alpha',p);

        % Display results
        if h_ai_acc1 == 0
            fprintf('ACC comparing anhedonia and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_ai_acc1);
        else
            fprintf('ACC comparing anhedonia and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_ai_acc1);
        end
        %% mOFC
        [p_gi_ofc1,h_gi_ofc1,stat_gi_ofc1] = ranksum(bootstat_gi_sim(:), bootstat_ofc_gi_rest,'Alpha',p);

        % Display results
        if h_gi_ofc1 == 0
            fprintf('mOFC comparing general distress and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_gi_ofc1);
        else
            fprintf('mOFC comparing general distress and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_gi_ofc1);
        end

        [p_ai_ofc1,h_ai_ofc1,stat_ai_ofc1] = ranksum(bootstat_ai_sim(:), bootstat_ofc_ai_rest,'Alpha',p);

        % Display results
        if h_ai_ofc1 == 0
            fprintf('mOFC comparing anhedonia and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_ai_ofc1);
        else
            fprintf('mOFC comparing anhedonia and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_ai_ofc1);
        end

        %% Amyg
        [p_gi_amyg1,h_gi_amyg1,stat_gi_amyg1] = ranksum(bootstat_gi_sim(:), bootstat_amyg_gi_rest,'Alpha',p);

        % Display results
        if h_gi_amyg1 == 0
            fprintf('Amygdala comparing general distress and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_gi_amyg1);
        else
            fprintf('Amygdala comparing general distress and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_gi_amyg1);
        end

        [p_ai_amyg1,h_ai_amyg1,stat_ai_amyg1] = ranksum(bootstat_ai_sim(:), bootstat_amyg_ai_rest,'Alpha',p);

        % Display results
        if h_ai_amyg1 == 0
            fprintf('Amygdala comparing anhedonia and inflammation is not significantly different from simulated values using resting state (p = %.4f).\n', p_ai_amyg1);
        else
            fprintf('Amygdala comparing anhedonia and inflammation is significantly different from simulated values using resting state (p = %.4f).\n', p_ai_amyg1);
        end

        %%% MID
        %% ACC
        [p_gi_acc,h_gi_acc,stat_gi_acc] = ranksum(bootstat_gi_sim(:), bootstat_acc_gi_mid,'Alpha',p);

        % Display results
        if h_gi_acc == 0
            fprintf('ACC comparing general distress and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_gi_acc);
        else
            fprintf('ACC comparing general distress and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_gi_acc);
        end

        [p_ai_acc,h_ai_acc,stat_ai_acc] = ranksum(bootstat_ai_sim(:), bootstat_acc_ai_mid,'Alpha',0.01/12);

        % Display results
        if h_ai_acc == 0
            fprintf('ACC comparing anhedonia and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_ai_acc);
        else
            fprintf('ACC comparing anhedonia and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_ai_acc);
        end
        %% mOFC
        [p_gi_ofc,h_gi_ofc,stat_gi_ofc] = ranksum(bootstat_gi_sim(:), bootstat_ofc_gi_mid,'Alpha',p);

        % Display results
        if h_gi_ofc == 0
            fprintf('mOFC comparing general distress and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_gi_ofc);
        else
            fprintf('mOFC comparing general distress and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_gi_ofc);
        end

        [p_ai_ofc,h_ai_ofc,stat_ai_ofc] = ranksum(bootstat_ai_sim(:), bootstat_ofc_ai_mid,'Alpha',p);

        % Display results
        if h_ai_ofc == 0
            fprintf('mOFC comparing anhedonia and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_ai_ofc);
        else
            fprintf('mOFC comparing anhedonia and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_ai_ofc);
        end

        %% Amyg
        [p_gi_amyg,h_gi_amyg,stat_gi_amyg] = ranksum(bootstat_gi_sim(:), bootstat_amyg_gi_mid,'Alpha',p);

        % Display results
        if h_gi_amyg == 0
            fprintf('Amygdala comparing general distress and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_gi_amyg);
        else
            fprintf('Amygdala comparing general distress and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_gi_amyg);
        end

        [p_ai_amyg,h_ai_amyg,stat_ai_amyg] = ranksum(bootstat_ai_sim(:), bootstat_amyg_ai_mid,'Alpha',p);

        % Display results
        if h_ai_amyg == 0
            fprintf('Amygdala comparing anhedonia and inflammation is not significantly different from simulated values using MID (p = %.4f).\n', p_ai_amyg);
        else
            fprintf('Amygdala comparing anhedonia and inflammation is significantly different from simulated values using MID (p = %.4f).\n', p_ai_amyg);
        end
        
        % plot results, only running this if I don't also want to plot
        % visual control results at the same time
        if control_visual == 0
            figure();
            subplot(3,1,1);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_amyg_gi_mid); 
            title('Similarity of General Distress and Inflammation: MID amygdala transformations')
            subplot(3,1,2);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_ofc_gi_mid); 
            title('Similarity of General Distress and Inflammation: MID mOFC transformations')
            subplot(3,1,3);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_acc_gi_mid);
            title('Similarity of General Distress and Inflammation: MID ACC transformations')
    
            figure();
            subplot(3,1,1);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_amyg_gi_rest); 
            title('Similarity of General Distress and Inflammation: Resting state amygdala transformations')
            subplot(3,1,2);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_ofc_gi_rest); 
            title('Similarity of General Distress and Inflammation: Resting state mOFC transformations')
            subplot(3,1,3);
            histogram(bootstat_gi_sim); hold on; histogram(bootstat_acc_gi_rest); 
            title('Similarity of General Distress and Inflammation: Resting state ACC transformations')
    
            figure();
            subplot(3,1,1);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_amyg_ai_mid); 
            title('Similarity of Anhedonia-Apprehension and Inflammation: MID amygdala transformations')
            subplot(3,1,2);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_ofc_ai_mid); 
            title('Similarity of Anhedonia-Apprehension and Inflammation: MID mOFC transformations')
            subplot(3,1,3);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_acc_ai_mid); 
            title('Similarity of Anhedonia-Apprehension and Inflammation: MID ACC transformations')
    
            figure();
            subplot(3,1,1);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_amyg_ai_rest); 
            title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state amygdala transformations')
            subplot(3,1,2);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_ofc_ai_rest);
            title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state mOFC transformations')
            subplot(3,1,3);
            histogram(bootstat_ai_sim); hold on; histogram(bootstat_acc_ai_rest);
            title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state ACC transformations')
        end
    end
    
    if control_visual == 1
        inflam = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_visual_control/inflammation/corr_transform.mat');
        gendis = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_visual_control/general_distress/corr_transform.mat');
        anhed = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_visual_control/anhedonia/corr_transform.mat');
        
        % separate out visual correlations with outcomes into voxelwise
        % estimates in single columns that match the size of the actual MID
        % cortico-striatal-amygdala corr maps
        for i = 1:14
            if i == 1
                inflam_corr = inflam.corr_transform{i}(:);
                anhed_corr = anhed.corr_transform{i}(:);
                gendis_corr = gendis.corr_transform{i}(:);
            else
                inflam_corr = [inflam_corr; inflam.corr_transform{i}(:)];
                anhed_corr = [anhed_corr; anhed.corr_transform{i}(:)];
                gendis_corr = [gendis_corr; gendis.corr_transform{i}(:)];
            end
        end

        % generate a matching 1000 permutations of correlations across
        % these multivariate patterns
        for perm = 1:100
            [bootstat_ai_vis(:,perm),~] = bootstrp(10,@corr,anhed_corr,inflam_corr);
            [bootstat_gi_vis(:,perm),~] = bootstrp(10,@corr,gendis_corr,inflam_corr);
            [bootstat_ag_vis(:,perm),~] = bootstrp(10,@corr,anhed_corr,gendis_corr);
        end
        
        % plot visual control with simulated control with real MID/rest
        figure();
        subplot(3,1,1);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_amyg_gi_mid); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: MID VS-amygdala transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')        
        subplot(3,1,2);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_ofc_gi_mid); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: MID VS-mOFC transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')
        subplot(3,1,3);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_acc_gi_mid); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: MID VS-ACC transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')

        figure();
        subplot(3,1,1);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_amyg_gi_rest); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: Resting state VS-amygdala transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')
        subplot(3,1,2);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_ofc_gi_rest); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: Resting state VS-mOFC transformations')
        legend('Simulated data', 'VS-mOFC data', 'Visual control')
        subplot(3,1,3);
        histogram(bootstat_gi_sim); hold on; histogram(bootstat_acc_gi_rest); hold on; histogram(bootstat_gi_vis);
        title('Similarity of General Distress and Inflammation: Resting state VS-ACC transformations')
        legend('Simulated data', 'VS-ACC data', 'Visual control')

        figure();
        subplot(3,1,1);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_amyg_ai_mid); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: MID VS-amygdala transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')
        subplot(3,1,2);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_ofc_ai_mid); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: MID VS-mOFC transformations')
        legend('Simulated data', 'VS-mOFC data', 'Visual control')
        subplot(3,1,3);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_acc_ai_mid); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: MID VS-ACC transformations')
        legend('Simulated data', 'VS-ACC data', 'Visual control')

        figure();
        subplot(3,1,1);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_amyg_ai_rest); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state VS-amygdala transformations')
        legend('Simulated data', 'VS-amygdala data', 'Visual control')
        subplot(3,1,2);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_ofc_ai_rest); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state VS-mOFC transformations')
        legend('Simulated data', 'VS-mOFC data', 'Visual control')
        subplot(3,1,3);
        histogram(bootstat_ai_sim); hold on; histogram(bootstat_acc_ai_rest); hold on; histogram(bootstat_ai_vis);
        title('Similarity of Anhedonia-Apprehension and Inflammation: Resting state VS-ACC transformations')
        legend('Simulated data', 'VS-ACC data', 'Visual control')

        %%% Resting state visual cortex comparison
        %% ACC
        [p_gi_acc1_vis,h_gi_acc1_vis,stat_gi_acc1_vis] = ranksum(bootstat_gi_vis(:), bootstat_acc_gi_rest,'Alpha',p);

        % Display results
        if h_gi_acc1_vis == 0
            fprintf('ACC comparing general distress and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_acc1_vis);
        else
            fprintf('ACC comparing general distress and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_acc1_vis);
        end

        [p_ai_acc1_vis,h_ai_acc1_vis,stat_ai_acc1_vis] = ranksum(bootstat_ai_vis(:), bootstat_acc_ai_rest,'Alpha',p);

        % Display results
        if h_ai_acc1_vis == 0
            fprintf('ACC comparing anhedonia and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_acc1_vis);
        else
            fprintf('ACC comparing anhedonia and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_acc1_vis);
        end
        %% mOFC
        [p_gi_ofc1_vis,h_gi_ofc1_vis,stat_gi_ofc1_vis] = ranksum(bootstat_gi_vis(:), bootstat_ofc_gi_rest,'Alpha',p);

        % Display results
        if h_gi_ofc1_vis == 0
            fprintf('mOFC comparing general distress and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_ofc1_vis);
        else
            fprintf('mOFC comparing general distress and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_ofc1_vis);
        end

        [p_ai_ofc1_vis,h_ai_ofc1_vis,stat_ai_ofc1_vis] = ranksum(bootstat_ai_vis(:), bootstat_ofc_ai_rest,'Alpha',0.01/12);

        % Display results
        if h_ai_ofc1_vis == 0
            fprintf('mOFC comparing anhedonia and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_ofc1_vis);
        else
            fprintf('mOFC comparing anhedonia and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_ofc1_vis);
        end

        %% Amyg
        [p_gi_amyg1_vis,h_gi_amyg1_vis,stat_gi_amyg1_vis] = ranksum(bootstat_gi_vis(:), bootstat_amyg_gi_rest,'Alpha',p);

        % Display results
        if h_gi_amyg1_vis == 0
            fprintf('Amygdala comparing general distress and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_amyg1_vis);
        else
            fprintf('Amygdala comparing general distress and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_gi_amyg1_vis);
        end

        [p_ai_amyg1_vis,h_ai_amyg1_vis,stat_ai_amyg1_vis] = ranksum(bootstat_ai_vis(:), bootstat_amyg_ai_rest,'Alpha',p);

        % Display results
        if h_ai_amyg1_vis == 0
            fprintf('Amygdala comparing anhedonia and inflammation is not significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_amyg1_vis);
        else
            fprintf('Amygdala comparing anhedonia and inflammation is significantly different from visual cortex values using resting state (p = %.4f).\n', p_ai_amyg1_vis);
        end

        %%% MID visual cortex comparison
        %% ACC
        [p_gi_acc_vis,h_gi_acc_vis,stat_gi_acc_vis] = ranksum(bootstat_gi_vis(:), bootstat_acc_gi_mid,'Alpha',p);

        % Display results
        if h_gi_acc_vis == 0
            fprintf('ACC comparing general distress and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_acc);
        else
            fprintf('ACC comparing general distress and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_acc);
        end

        [p_ai_acc_vis,h_ai_acc_vis,stat_ai_acc_vis] = ranksum(bootstat_ai_vis(:), bootstat_acc_ai_mid,'Alpha',0.01/12);

        % Display results
        if h_ai_acc_vis == 0
            fprintf('ACC comparing anhedonia and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_acc_vis);
        else
            fprintf('ACC comparing anhedonia and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_acc_vis);
        end
        %% mOFC
        [p_gi_ofc_vis,h_gi_ofc_vis,stat_gi_ofc_vis] = ranksum(bootstat_gi_vis(:), bootstat_ofc_gi_mid,'Alpha',p);

        % Display results
        if h_gi_ofc_vis == 0
            fprintf('mOFC comparing general distress and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_ofc_vis);
        else
            fprintf('mOFC comparing general distress and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_ofc_vis);
        end

        [p_ai_ofc_vis,h_ai_ofc_vis,stat_ai_ofc_vis] = ranksum(bootstat_ai_vis(:), bootstat_ofc_ai_mid,'Alpha',p);

        % Display results
        if h_ai_ofc_vis == 0
            fprintf('mOFC comparing anhedonia and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_ofc_vis);
        else
            fprintf('mOFC comparing anhedonia and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_ofc_vis);
        end

        %% Amyg
        [p_gi_amyg_vis,h_gi_amyg_vis,stat_gi_amyg_vis] = ranksum(bootstat_gi_vis(:), bootstat_amyg_gi_mid,'Alpha',p);

        % Display results
        if h_gi_amyg_vis == 0
            fprintf('Amygdala comparing general distress and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_amyg_vis);
        else
            fprintf('Amygdala comparing general distress and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_gi_amyg_vis);
        end

        [p_ai_amyg_vis,h_ai_amyg_vis,stat_ai_amyg_vis] = ranksum(bootstat_ai_vis(:), bootstat_amyg_ai_mid,'Alpha',p);

        % Display results
        if h_ai_amyg_vis == 0
            fprintf('Amygdala comparing anhedonia and inflammation is not significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_amyg_vis);
        else
            fprintf('Amygdala comparing anhedonia and inflammation is significantly different from visual cortex values using MID (p = %.4f).\n', p_ai_amyg_vis);
        end
        
    end
    
end

if display_regressors == 1
    mid = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/regressors_mid.mat');
    rest = load('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/regressors_rest.mat');

    stats1 = grpstats(mid.regressors, [], {'mean', 'median', 'std', 'min', 'max', 'iqr'});
    stats2 = grpstats(rest.regressors, [], {'mean', 'median', 'std', 'min', 'max', 'iqr'});
    
    display(stats1)
    display(stats2)

    figure();
    %subplot(2,1,1)
    h1 = heatmap(corr([mid.regressors.inflammation,mid.regressors.longGeneralDistress,mid.regressors.longAnhedonia]));
    title('Correlations of clinical outcomes: MID task sample')
    h1.XDisplayLabels = {'Inflammation','long. General Distress','long. Anhedonia-Apprehension'};
    h1.YDisplayLabels = {'Inflammation','long. General Distress','long. Anhedonia-Apprehension'};
    %subplot(2,1,2)
    figure();
    h2 = heatmap(corr([rest.regressors.inflammation,rest.regressors.longGeneralDistress,rest.regressors.longAnhedonia]));
    title('Correlations of clinical outcomes: Resting state sample')
    h2.XDisplayLabels = {'Inflammation','long. General Distress','long. Anhedonia-Apprehension'};
    %h2.YDisplayLabels = {' ',' ',' '};
    h2.YDisplayLabels = {'Inflammation','long. General Distress','long. Anhedonia-Apprehension'};
    
    mdl1 = fitlm(mid.regressors,'inflammation~longGeneralDistress');
    mdl2 = fitlm(mid.regressors,'inflammation~longAnhedonia');

    figure();
    subplot(1,2,1);
    scatter(mid.regressors.inflammation, mid.regressors.longGeneralDistress);
    title('Inflammation and General Distress: MID task sample')
    subplot(1,2,2);
    scatter(mid.regressors.inflammation, mid.regressors.longAnhedonia);
    title('Inflammation and Anhedonia-Apprehension: MID task sample')

    figure();
    subplot(1,2,1);
    scatter(mid.regressors.inflammation, mid.regressors.longGeneralDistress);
    title('Inflammation and General Distress: Resting state sample')
    subplot(1,2,2);
    scatter(mid.regressors.inflammation, mid.regressors.longAnhedonia);
    title('Inflammation and Anhedonia-Apprehension: Resting state sample')

end
