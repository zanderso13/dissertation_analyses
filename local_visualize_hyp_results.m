
mid = 0;
rest = 1;
visualize_pls_results = 1;

if visualize_pls_results == 1
    if mid == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/mid_shortcuts/inflammation/')
    end
    if rest == 1
        cd('/Users/zacharyanderson/Documents/ACNlab/dissertation_final/rest_shortcuts/inflammation')
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

        loadings_to_visualize(:,:,i) = reshape(transformXL{i}(:,1),[828,71]);
    end

    % visualize MSE
    % no differences in MSE. Distributions look super similar
    figure(); histogram(allMSE(:,1)); hold on; histogram(allMSE(:,3)); hold on; xline(allMSE(1,2),'LineWidth', 5); hold on; histogram(allMSE(:,4));
    figure(); histogram(allMSE(:,1)); hold on; xline(allMSE(1,2),'LineWidth', 5); hold on; histogram(allMSE(:,4));
    p_aligned_pls_mse = sum(allMSE(:,1)>allMSE(1,2))./100; % 
    p_transforms_pls_mse = sum(allMSE(:,4)>allMSE(1,2))./100;

    % visualize pct variance covered
    % unaligned data is no better than the null dist. Aligned data is WAY
    % better than both in terms of percent variance accounted for
    pctall = figure(); histogram(allPCT(:,1)); hold on; histogram(allPCT(:,3)); hold on; xline(allPCT(1,2),'LineWidth', 5); hold on; histogram(allPCT(:,4));
    pcthyp = figure(); histogram(allPCT(:,1)); hold on; xline(allPCT(1,2),'LineWidth', 5); hold on; histogram(allPCT(:,4));
    savefig(pctall,'pctall.fig'); savefig(pcthyp,'pct_hyperalignment.fig');
    p_aligned_pls_pct = sum(allPCT(:,1)<allPCT(1,2))./100; % p < 0.01
    p_transform_pls_pct = sum(allPCT(:,4)<allPCT(1,2))./100; % p < 0.01
    
    brainfile = filenames('hypatl*mat');
    load(brainfile{1})

    figure(); montage(hypatl)
     
    

end