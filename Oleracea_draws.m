% Life cycle parameters draws per population
function draws = Oleracea_draws(run, population,populations,Histo_sizes,...
    sensi_para,draws,population_year,Do_sensitivity,max_diameter,sensi,histo_length)

%% Stochasticity among years
PSB_mean = populations.(genvarname(['population','_',int2str(population)])).PSB_mean; 
Gs_mean_undis =populations.(genvarname(['population','_',int2str(population)])).Gs_mean_all(1,:);
Gs_mean_dis =populations.(genvarname(['population','_',int2str(population)])).Gs_mean_all(2,:);
Ss_mean_undis =populations.(genvarname(['population','_',int2str(population)])).Ss_mean_all(1,:);
Ss_mean_dis =populations.(genvarname(['population','_',int2str(population)])).Ss_mean_all(2,:);
Ss_Dis_std = populations.(genvarname(['population','_',int2str(population)])).Ss_std_all(2,:);
Adult_surv_mean =  populations.(genvarname(['population','_',int2str(population)])).Adult_surv_mean_all;
Adult_increase_mean = populations.(genvarname(['population','_',int2str(population)])).Adult_increase_mean_all;
Flowering_like_mean = populations.(genvarname(['population','_',int2str(population)])).Flowering_like_mean_all;
Fecundity_mean_fn = populations.(genvarname(['population','_',int2str(population)])).Fecundity_mean_fn_all;
Fecundity_mean_c = populations.(genvarname(['population','_',int2str(population)])).Fecundity_mean_c_all;

%% Seedbank
% function needs to be within limits 0-100%
SBincorp = PSB_mean(1);
if SBincorp < 0
    SBincorp = 0;
end
if SBincorp > 1
    SBincorp = 1;
end
% function needs to be within limits 0-100%
SBrate = PSB_mean(2);
if SBrate < 0
    SBrate = 0;
end
if SBrate > 1
    SBrate = 1;
end
clear PSB_mean
clear PSB_std
%% Emergence and seedling survival
for f= 1:1:length(Gs_mean_dis)
    % function needs to be within limits 0-100%
    Gs_dis(f) = Gs_mean_dis(f);
    if Gs_dis(f) < 0
        Gs_dis(f) = 0;
    end
    if Gs_dis(f) > 1
        Gs_dis(f) = 1;
    end
end
clear GS_mean_dis
clear GS_std_dis
for f= 1:1:length(Gs_mean_undis)
    % function needs to be within limits 0-100%
    Gs_undis(f) = Gs_mean_undis(f);
    if Gs_undis(f) < 0
        Gs_undis(f) = 0;
    end
    if Gs_undis(f) > 1
        Gs_undis(f) = 1; %#ok<*AGROW>
    end
end
clear addition
clear GS_mean_dis
clear GS_std_dis
%addition = randn;
for f= 1:1:length(Ss_mean_dis)
    % function needs to be within limits 0-100%
    Ss_dis(f) = Ss_mean_dis(f);%+ (Ss_Dis_std(f).*addition);
    if Ss_dis(f) < 0
        Ss_dis(f) = 0;
    end
    if Ss_dis(f) > 1
        Ss_dis(f) = 1;
    end
end
clear Ss_mean_dis
clear Ss_std_dis
Ss_dis(5) = Ss_dis(3);
for f= 1:1:length(Ss_mean_undis)
    % function needs to be within limits 0-100%
    Ss_undis(f) = Ss_mean_undis(f);
    if Ss_undis(f) < 0
        Ss_undis(f) = 0;
    end
    if Ss_undis(f) > 1
        Ss_undis(f) = 1;
    end
end
clear addition
clear Ss_mean_dis
clear Ss_std_dis
Ss_undis(5) = Ss_undis(3);

%% Adult survival likelihoods
s_max = 2;
% function needs to be within limit below 100%
while s_max > 1
    Bc = Adult_surv_mean(population_year,1);
    B1 = Adult_surv_mean(population_year,2);
    B2 = Adult_surv_mean(population_year,3);
    for i = 1:1:max_diameter
        x = log10(Histo_sizes(i,1)+1);
        Adult_surv(2,i) = Bc + (B1*(x)) +  (B2*(x^2));
        Adult_surv(1,i) = Histo_sizes(i,1);
        
        if  Adult_surv(2,i) < 0
            Adult_surv(2,i) = 0;
        end
        if  Adult_surv(2,i) >1
            Adult_surv(2,i) = 1;
        end
    end
    s_max = max(Adult_surv(2,:));
end
clear All not needed % shortened to save text
%% Adult size increase
s_max = 2; % think about limits, if applicable
while s_max >1
    Bc = Adult_increase_mean(population_year,1);
    B1 =Adult_increase_mean(population_year,2);
    for i = 1:1:max_diameter
        x = log10((Histo_sizes(i,1))+1);
        Adult_increase(2,i) = (Bc + (B1*(x)))+1; % Note that % increase rate is transforned to % rate (so + 1)
        Adult_increase(1,i) = Histo_sizes(i,1);
    end
    s_max = 0.5;
end
clear All not needed % shortened to save text
%% Likelihood of flowering
% function needs to be within limits of 0 and 100%
s_max = 2;
s_min = -1;
while s_max > 1 || s_min <0
    Bc = Flowering_like_mean(population_year,1);
    B1 = Flowering_like_mean(population_year,2);
    for i = 1:1:max_diameter
        x = Histo_sizes(i,1);
        Flowering_like(2,i) = 1/(1+(exp(Bc+(B1.*x))));
        Flowering_like(1,i) = x;
    end
    s_max = max(Flowering_like(2,:));
    s_min = min(Flowering_like(2,:));
end
clear All not needed % shortened to save text
%% Fecundity
% Function has a lower limit of 0 and no maximum
s_min = -1;
while s_min <0
    for i = 1:1:max_diameter
        x = Histo_sizes(i,1);
        Fecundity(1,i) = Histo_sizes(i,1);
        % Constant function
        Fecundity(2,i) = -1;
        while Fecundity(2,i) < 0
            Fecundity(2,i) = Fecundity_mean_c(population_year);
        end       
    end
    s_min = min(Fecundity(2,:));
end

%% Sensitivity
if exist('Do_sensitivity') ~= 0
    if strcmpi('y',Do_sensitivity) ~= 0
        Oleracea_Sensitivy_runss;
    end
end

%% Make it arrays per population
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Flowering_like = Flowering_like;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Fecundity = Fecundity;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).SBincorp = SBincorp;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).SBrate = SBrate;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Gs_undis = Gs_undis;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Ss_undis = Ss_undis;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Gs_dis = Gs_dis;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Ss_dis = Ss_dis;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Adult_surv = Adult_surv;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Adult_increase = Adult_increase;
draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).histo_length = histo_length;



