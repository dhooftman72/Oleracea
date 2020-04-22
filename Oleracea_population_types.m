function  population_types = Oleracea_population_types
load('parameters.mat')
%load('parameters_start.mat')
for population_type = 1:1:3
    if population_type ==1
        population_types.Kim.PSB_mean = PSB (:,1);
        population_types.Kim.PSB_std =   PSB (:,2);
        
        population_types.Kim.Gs_mean_all(1,:) = Gs_Undis_mean;
        population_types.Kim.Gs_std_all(1,:) = Gs_Undis_std;
        population_types.Kim.Ss_mean_all(1,:) = Ss_Undis_mean;
        population_types.Kim.Ss_std_all(1,:) = Ss_Undis_std;
        population_types.Kim.Gs_mean_all(2,:) = Gs_Dis_mean;
        population_types.Kim.Gs_std_all(2,:) = Gs_Dis_std;
        population_types.Kim.Ss_mean_all(2,:) = Ss_Dis_mean;
        population_types.Kim.Ss_std_all(2,:) = Ss_Dis_std;
        
        population_types.Kim.histo_org_length= (size(Histo_Base_K,1));
        population_types.Kim.Histo_base(1:(population_types.Kim.histo_org_length(1))) = Histo_Base_K(:,2);
        population_types.Kim.Adult_surv_mean_all = Adult_surv_mean_K(1:9,:);
        population_types.Kim.Adult_surv_std_all = Adult_surv_std_K(1:9,:);
        population_types.Kim.Adult_increase_mean_all = Adult_increase_mean_K(1:9,:);
        population_types.Kim.Adult_increase_std_all =Adult_increase_std_K(1:9,:);
        population_types.Kim.Flowering_like_mean_all = Flowering_like_mean_K(1:9,:);
        population_types.Kim.Flowering_like_std_all = Flowering_like_std_K(1:9,:);
        population_types.Kim.Fecundity_mean_fn_all = Fecundity_mean_K(1:7,1);
        population_types.Kim.Fecundity_std_fn_all = Fecundity_std_K(1:7,1);
        population_types.Kim.Fecundity_mean_c_all = Fecundity_mean_K(1:7,2);
        population_types.Kim.Fecundity_std_c_all = Fecundity_std_K(1:7,2);
        
        % We include no std
        population_types.Kim.PSB_std(:,:) =  0;
        population_types.Kim.Gs_std_all(:,:) =  0;
        population_types.Kim.Ss_std_all(:,:) =  0;
        population_types.Kim.Adult_surv_std_all(:,:) =  0;
        population_types.Kim.Adult_increase_std_all(:,:) =  0;
        population_types.Kim.Flowering_like_std_all(:,:) =  0;
        population_types.Kim.Fecundity_std_fn_all(:,:) =  0;
        population_types.Kim.Fecundity_std_c_all(:,:) =  0;

    elseif population_type ==2
        population_types.OH.PSB_mean = PSB (:,1);
        population_types.OH.PSB_std =   PSB (:,2);
        
        population_types.OH.Gs_mean_all(1,:) = Gs_Undis_mean;
        population_types.OH.Gs_std_all(1,:) = Gs_Undis_std;
        population_types.OH.Ss_mean_all(1,:) = Ss_Undis_mean;
        population_types.OH.Ss_std_all(1,:) = Ss_Undis_std;
        population_types.OH.Gs_mean_all(2,:) = Gs_Dis_mean;
        population_types.OH.Gs_std_all(2,:) = Gs_Dis_std;
        population_types.OH.Ss_mean_all(2,:) = Ss_Dis_mean;
        population_types.OH.Ss_std_all(2,:) = Ss_Dis_std;
        
        population_types.OH.histo_org_length= (size(Histo_Base_OH,1));
        population_types.OH.Histo_base(1:(population_types.OH.histo_org_length(1))) = Histo_Base_OH(:,2);
        population_types.OH.Adult_surv_mean_all = Adult_surv_mean_OH(1:9,:);
        population_types.OH.Adult_surv_std_all = Adult_surv_std_OH(1:9,:);
        population_types.OH.Adult_increase_mean_all = Adult_increase_mean_OH(1:9,:);
        population_types.OH.Adult_increase_std_all =Adult_increase_std_OH(1:9,:);
        population_types.OH.Flowering_like_mean_all = Flowering_like_mean_OH(1:9,:);
        population_types.OH.Flowering_like_std_all = Flowering_like_std_OH(1:9,:);
        population_types.OH.Fecundity_mean_fn_all = Fecundity_mean_OH(1:7,1);
        population_types.OH.Fecundity_std_fn_all = Fecundity_std_OH(1:7,1);
        population_types.OH.Fecundity_mean_c_all = Fecundity_mean_OH(1:7,2);
        population_types.OH.Fecundity_std_c_all = Fecundity_std_OH(1:7,2);
        
        % We include no std
            population_types.OH.PSB_std(:,:) =  0;
            population_types.OH.Gs_std_all(:,:) =  0;
            population_types.OH.Ss_std_all(:,:) =  0;
            population_types.OH.Adult_surv_std_all(:,:) =  0;
            population_types.OH.Adult_increase_std_all(:,:) =  0;
            population_types.OH.Flowering_like_std_all(:,:) =  0;
            population_types.OH.Fecundity_std_fn_all(:,:) =  0;
            population_types.OH.Fecundity_std_c_all(:,:) =  0;

        
    elseif population_type ==3
        population_types.Win.PSB_mean = PSB (:,1);
        population_types.Win.PSB_std =   PSB (:,2);
        
        population_types.Win.Gs_mean_all(1,:) = Gs_Undis_mean;
        population_types.Win.Gs_std_all(1,:) = Gs_Undis_std;
        population_types.Win.Ss_mean_all(1,:) = Ss_Undis_mean;
        population_types.Win.Ss_std_all(1,:) = Ss_Undis_std;
        population_types.Win.Gs_mean_all(2,:) = Gs_Dis_mean;
        population_types.Win.Gs_std_all(2,:) = Gs_Dis_std;
        population_types.Win.Ss_mean_all(2,:) = Ss_Dis_mean;
        population_types.Win.Ss_std_all(2,:) = Ss_Dis_std;
        
        population_types.Win.histo_org_length= (size(Histo_Base_W,1));
        population_types.Win.Histo_base(1:(population_types.Win.histo_org_length(1))) = Histo_Base_W(:,2);
        population_types.Win.Adult_surv_mean_all = Adult_surv_mean_W(1:9,:);
        population_types.Win.Adult_surv_std_all = Adult_surv_std_W(1:9,:);
        population_types.Win.Adult_increase_mean_all = Adult_increase_mean_W(1:9,:);
        population_types.Win.Adult_increase_std_all =Adult_increase_std_W(1:9,:);
        population_types.Win.Flowering_like_mean_all = Flowering_like_mean_W(1:9,:);
        population_types.Win.Flowering_like_std_all = Flowering_like_std_W(1:9,:);
        population_types.Win.Fecundity_mean_fn_all = Fecundity_mean_W(1:7,1);
        population_types.Win.Fecundity_std_fn_all = Fecundity_std_W(1:7,1);
        population_types.Win.Fecundity_mean_c_all = Fecundity_mean_W(1:7,2);
        population_types.Win.Fecundity_std_c_all = Fecundity_std_W(1:7,2);
        
        % We include no std
            population_types.Win.PSB_std(:,:) =  0;
            population_types.Win.Gs_std_all(:,:) =  0;
            population_types.Win.Ss_std_all(:,:) =  0;
            population_types.Win.Adult_surv_std_all(:,:) =  0;
            population_types.Win.Adult_increase_std_all(:,:) =  0;
            population_types.Win.Flowering_like_std_all(:,:) =  0;
            population_types.Win.Fecundity_std_fn_all(:,:) =  0;
            population_types.Win.Fecundity_std_c_all(:,:) =  0;
    end
end