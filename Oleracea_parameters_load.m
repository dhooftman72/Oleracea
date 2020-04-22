function Oleracea_parameters_load
%% Paraqmeters: get all the functions we want first
    cd Data_analyses_results
    load('seedbank_functions_results.mat');
    load('sowing_exp_functions_results.mat');
    load('survival_functions_results.mat');
    load('size_functions_results.mat');
    load('flowering_functions_results.mat');
    load('fecundity_functions_results.mat');
    load('starting_histograms.mat')
    cd ../
    
    % Seedbank
    PSB(1,1) = Seedbank_function(10);
    PSB(1,2) = Seedbank_function(11);
    PSB(2,1) = Seedbank_function(16);
    PSB(2,2) = Seedbank_function(17);
    
    % Germination, Can be adapted question dependent!!!
        % Disturbed
    temp = Sowing_Patch_Remain([1:24,49:57],:);
    temp(:,6) = Sowing_Germination_Patch([1:24,49:57],6);
    for i = 1:1:6
       temp2 = temp(:,i);
       temp2(isnan(temp2) == 1) = [];
       Gs_Dis_mean(i) = mean(temp2);
       Gs_Dis_std(i) = std(temp2);
       clear temp2
    end
    
    clear temp
        % disturbed
        temp = Sowing_Patch_Remain([25:48,58:66],:);
         temp(:,6) = Sowing_Germination_Patch([25:48,58:66],6);
    for i = 1:1:6
       temp2 = temp(:,i);
       temp2(isnan(temp2) == 1) = [];
       Gs_Undis_mean(i) = mean(temp2);
       Gs_Undis_std(i) = std(temp2);
       clear temp2
    end
    clear temp
    
    %% Seedling Survival  Can be adapted question dependent!!!
        %Disturbed
    temp = Sowing_survival_Patch([1:24,49:57],[1:4]);
    temp(:,5) = Sowing_survival_Patch([1:24,49:57],6);
    for i = 1:1:5
       temp2 = temp(:,i);
       temp2(isnan(temp2) == 1) = [];
       Ss_Dis_mean(i) = mean(temp2);
       Ss_Dis_std(i) = std(temp2);
       clear temp2
    end
    clear temp
        % Undisturbed
        temp = Sowing_survival_Patch([25:48,58:66],[1:4]);
         temp(:,5) = Sowing_Germination_Patch([25:48,58:66],6);
    for i = 1:1:5
       temp2 = temp(:,i);
       temp2(isnan(temp2) == 1) = [];
       Ss_Undis_mean(i) = mean(temp2);
       Ss_Undis_std(i) = std(temp2);
       clear temp2
    end
    clear temp
    
    % Likelihood of adults survival as 2nd order polymial log function
    Adult_surv_mean_K = Survival_quad_function_K(:,[5:7]);
    Adult_surv_std_K = Survival_quad_function_K(:,[8:10]);
    Adult_surv_mean_OH = Survival_quad_function_OH(:,[5:7]);
    Adult_surv_std_OH = Survival_quad_function_OH(:,[8:10]);
    Adult_surv_mean_W = Survival_quad_function_W(:,[5:7]);
    Adult_surv_std_W = Survival_quad_function_W(:,[8:10]);
    
    % Size increae as logaritmic function
    Adult_increase_mean_K = Growth_function_K(:,[4,5]); % Note that this is increase only, so 0 is equal size, current size needs to be added.
    Adult_increase_std_K = Growth_function_K(:,[6,7]);
    Adult_increase_mean_OH = Growth_function_OH(:,[4,5]); % Note that this is increase only, so 0 is equal size, current size needs to be added.
    Adult_increase_std_OH = Growth_function_OH(:,[6,7]);
    Adult_increase_mean_W = Growth_function_W(:,[4,5]); % Note that this is increase only, so 0 is equal size, current size needs to be added.
    Adult_increase_std_W = Growth_function_W(:,[6,7]);
      
    % Likelihood flowering as Logit function
    Flowering_like_mean_K = Flowering_function_K(:,[1:2]);
    Flowering_like_std_K = Flowering_function_K(:,[5:6]);
    Flowering_like_mean_OH = Flowering_function_OH(:,[1:2]);
    Flowering_like_std_OH = Flowering_function_OH(:,[5:6]);
    Flowering_like_mean_W = Flowering_function_W(:,[1:2]);
    Flowering_like_std_W = Flowering_function_W(:,[5:6]);
    
    % Fecundity as linear function with 0,0 intercept
    Fecundity_mean_K = Seeds_no_function_K(:,4);
    Fecundity_std_K = Seeds_no_function_K(:,6);
    Fecundity_mean_OH = Seeds_no_function_OH(:,4);
    Fecundity_std_OH = Seeds_no_function_OH(:,6);
    Fecundity_mean_W = Seeds_no_function_W(:,4);
    Fecundity_std_W = Seeds_no_function_W(:,6);
    
    % Fecundity as no (constant) function 
    Fecundity_mean_K(:,2) =  Seeds_no_function_K(:,13);
    Fecundity_std_K(:,2) =  Seeds_no_function_K(:,10);    
    Fecundity_mean_OH(:,2) =  Seeds_no_function_OH(:,13); 
    Fecundity_std_OH(:,2) =  Seeds_no_function_OH(:,10);  
    Fecundity_mean_W(:,2) =  Seeds_no_function_W(:,13);
    Fecundity_std_W(:,2) =  Seeds_no_function_W(:,10); 
    
    save('parameters', 'PSB', 'Gs_Dis_mean', 'Gs_Dis_std','Gs_Undis_mean','Gs_Undis_std',...
        'Ss_Dis_mean', 'Ss_Dis_std', 'Ss_Undis_mean', 'Ss_Undis_std', 'Adult_surv_mean_K', 'Adult_surv_std_K',...
        'Adult_surv_mean_OH', 'Adult_surv_std_OH', 'Adult_surv_mean_W', 'Adult_surv_std_W',...
        'Adult_increase_mean_K', 'Adult_increase_std_K', 'Adult_increase_mean_OH', 'Adult_increase_std_OH',...
        'Adult_increase_mean_W', 'Adult_increase_std_W', 'Flowering_like_mean_K', 'Flowering_like_std_K',...
        'Flowering_like_mean_OH', 'Flowering_like_std_OH', 'Flowering_like_mean_W', 'Flowering_like_std_W',...
        'Fecundity_mean_K', 'Fecundity_std_K', 'Fecundity_mean_OH', 'Fecundity_std_OH', 'Fecundity_mean_W',...
        'Fecundity_std_W','Histo_Base_K', 'Histo_Base_OH','Histo_Base_W');

