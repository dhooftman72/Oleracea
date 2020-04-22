function [run_max] =...
    Oleracea_actual_simulations(run,population_types,max_density,Prop_among,Sensis,...
    max_diameter,size_steps,sensi,histo_length, run_max,extention)
% Note: populations, Histo_sizes, Prop_among and Distance matrix are needed
% for embedded functions
rng('shuffle')
load('Prop_among.mat')
define_function

%% redefine histogram with a seedbank and 2 seedlings classes
Histo_base_model(1:3,pop_max) = 0;
Histo_base_model(4:histo_length,pop_max)= 0;
%% Initiate output parameters
for population = 1:1:pop_max
    lengths = populations.(genvarname(['population','_',int2str(population)])).histo_org_length;
    Histo_base_model(4:(3+lengths),population) = populations.(genvarname(['population','_',int2str(population)])).Histo_base;
    Start_individuals(population) = sum(Histo_base_model(4:histo_length,population));
    for i = 2:2:(histo_length-3)
        if Start_individuals(population) >0
            Existing_Stable_state((i/2),population) = (sum(Histo_base_model((i+2):(i+3),population)))./Start_individuals(population);
        else
            Existing_Stable_state((i/2),population) = 0;
        end
    end
end
% Base model Cabbage model per population per year
%% THIS IS THE REAL MODEL per year!!!!!!!!!!!!!!
% Make all random matrices first
adult_test = rand(year_max,pop_max);
seed_test = rand(year_max,pop_max);
class1_test = rand(year_max,pop_max);
class2_test= rand(year_max,pop_max);
reached = zeros(pop_max,1);
lambda_meta = 0;
occup_meta = 0;
for year = 1:1:year_max
    %display
    if year/sensi_para.show == round(year/sensi_para.show)
        clc
        display(sensi)
        display(run)
        display(year)
        if exist('Do_sensitivity') ~= 0
            if strcmpi('y',Do_sensitivity) ~= 0
                display('Sensitivity_mode')
            end
        end
    end
    % Pick the year of the parameters
    if strcmpi('y',Do_sensitivity) ~= 0
        population_year = population_years(year);
    else
        if strcmpi('y',sensi_para.all_years) ~= 0
            population_year = 7;
        else
            population_year = 0;
            while population_year == 0
                population_year = randi(6);
            end
        end
    end  
    %% From here per population!!!
    %% Initiate parameters and variables
    draws.org = {[1,2]};
    Colonisation = zeros(pop_max,1);
    Extinction = zeros(pop_max,1);
    for population = 1:1:pop_max
        %Initiate variables
        if year == 1
            Histo_individuals(:,population) = Histo_base_model(:,population)'; %#ok<*AGROW>
            Individuals(1,population) =  Start_individuals(population);
            total_individual(1) = sum(Individuals);
            Individuals_all = sum(round(Histo_individuals(2:histo_length,population)));
            Occupancy_adult(1) = (length(find(Start_individuals> sensi_para.threshold)))./pop_max;
            Occupancy_seeds(1) = (length(find( Histo_individuals(1) > sensi_para.threshold)))./pop_max;
            if Individuals_all > 0 % existing populations
                mortality(population) = 0;
                lambda(population) = 1;
                Colonisation(population) = 1;
            elseif Individuals_all == 0 % Non existing populations
                mortality(population) = 1;
                lambda(population) = 0;
            end
        end 
        % draw individual function parameters
        % here the function draws are made per year and sensitivities within
        draws = Oleracea_draws(run, population,populations,Histo_sizes,...
            sensi_para,draws,population_year,Do_sensitivity,max_diameter,sensi,histo_length);
    end
    
    %% Fecundity part (autumn)
    for population = 1:pop_max
        Fecundity = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Fecundity;
        Flowering_like = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Flowering_like;
        Flowering_plants = Histo_individuals(4:histo_length,population).* Flowering_like(2,:)'; %year t-1
        Seeds_total(population) = sum(Flowering_plants.*Fecundity(2,:)'); %year t-1
        % Round seeds to Threshold numbers: Make sure there are not infitive small number of seeds produced
        Seeds_total(population)= round(Seeds_total(population).*(1/sensi_para.threshold)).*sensi_para.threshold;
    end
    
    %% add seed flow from outside minus seed flow to outside
    %Oleracea_Seed_flow_external;
    if  strcmpi('y',sensi_para.wind_dispersal_present) ~= 0 || strcmpi('y',sensi_para.human_dispersal_present) ~= 0
        Oleracea_Seed_flow_year;
    end
    %% Growth part (spring)
    for population = 1:pop_max
        Histo_stop = 0;
        SBincorp = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).SBincorp;
        SBrate = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).SBrate;
        Gs_undis = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Gs_undis;
        Ss_undis = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Ss_undis;
        Gs_dis = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Gs_dis;
        Ss_dis = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Ss_dis;
        Adult_surv = (draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Adult_surv).*(1-sensi_para.mortality_disturbance_Adult);
        Adult_increase = draws.(genvarname(['r',int2str(run)])).(genvarname(['p',int2str(population)])).Adult_increase;
        Histo_individuals_meta = Histo_individuals(:,population);
        Seed_total = Seeds_total(population).*SBincorp; %year t-1
        % Seedbank survival from t-1 to t
        seedbank = Histo_individuals_meta(1).* SBrate;
        % Split between Disturbed and Undisturbed
        Histo_tmp_dis = Histo_individuals_meta(1:3).* sensi_para.proportion_disturbance;
        Histo_tmp_undis = Histo_individuals_meta(1:3).* (1- sensi_para.proportion_disturbance);
        Seed_total_dis = Seed_total.* sensi_para.proportion_disturbance;
        Seed_total_undis = Seed_total.* (1- sensi_para.proportion_disturbance);
        
        %Direct germination and to seedbank, note that the seeddling
        %survival is only present in the first year seedling component.
        % Germination at t
        % second year seedlings coming directly from first year in t-1
        % second year seedlings moving in the adult loop
        %DISTURBED
        to_first_year_dis = ((Gs_dis(1).*Ss_dis(1)).* Seed_total_dis)+ (Gs_dis(2).*Seed_total_dis); %
        to_seedbank_dis = Seed_total_dis- ((Gs_dis(1).* Seed_total_dis)+ (Gs_dis(2).*Seed_total_dis));
        germ_dis = (((Gs_dis(3)+Gs_dis(5))./2)+Gs_dis(4)).*Histo_tmp_dis(1);
        to_second_year_dis = Histo_tmp_dis(2).*(Ss_dis(2).*Ss_dis(3));
        to_adult_loop_dis = Histo_tmp_dis(3).*(Ss_dis(4).*Ss_dis(5));
        %UNDISTURBED
        to_first_year_undis = ((Gs_undis(1).*Ss_undis(1)).* Seed_total_undis)+ (Gs_undis(2).*Seed_total_undis); %
        to_seedbank_undis = Seed_total_undis- ((Gs_undis(1).* Seed_total_undis)+ (Gs_undis(2).*Seed_total_undis));
        germ_undis = (((Gs_undis(3)+Gs_undis(5))./2)+Gs_undis(4)).*Histo_tmp_undis(1);
        to_second_year_undis = Histo_tmp_undis(2).*(Ss_undis(2).*Ss_undis(3));
        to_adult_loop_undis = Histo_tmp_undis(3).*(Ss_undis(4).*Ss_undis(5));
        %JOIN
        to_first_year = to_first_year_dis + to_first_year_undis;
        to_seedbank = to_seedbank_dis +  to_seedbank_undis;
        germ =  germ_dis +  germ_undis;
        to_second_year = to_second_year_dis + to_second_year_undis;
        to_adult_loop = to_adult_loop_dis + to_adult_loop_undis;
        clear Histo_tmp_dis
        clear Histo_tmp_undis
        clear Seed_total_dis
        clear Seed_total_undis
        clear to_first_year_dis
        clear to_first_year_undis
        clear to_seedbank_dis
        clear to_seedbank_undis;
        clear germ_dis
        clear germ_undis;
        clear to_second_year_dis
        clear to_second_year_undis;
        clear to_adult_loop_dis
        clear to_adult_loop_undis;
        
        %subtract Germ from seed bank and add new seeds
        seedbank = seedbank-germ; % is seedbank at t
        seedbank = seedbank + to_seedbank;
        
        % add germination to first year seedlings
        to_first_year = to_first_year + germ; % first year seedlings at t
        
        %within adult_loop
        Adult_plants = zeros(4,(histo_length-3));
        Adult_plants(1,:) = Adult_increase(1,:);
        Adult_plants(2,:) = Histo_individuals_meta(4:histo_length).*  (Adult_surv(2,:)'); %year t-1
        Adult_plants(3,:) =  floor(Adult_plants(1,:).*Adult_increase(2,:));
        Adult_plants(4,:) = 0;
        
        % restrictions
        Adult_plants(3, Adult_plants(3,:)>100) = max_diameter;
        Adult_plants(3, Adult_plants(3,:)<2) = 2;
        
        for f = 1:1:max_diameter
            Adult_plants(4, (Adult_plants(3,f)-1)) = Adult_plants(2,f) + Adult_plants(4, (Adult_plants(3,f)-1));
        end
        Adult_plants(2,:) =  Adult_plants(4,:);
        Adult_plants(4,:) = [];
        Adult_plants(3,:) = [];
        
        
        % collate all transitions
        Histo_individuals_meta(4:histo_length) = Adult_plants(2,:);
        Histo_individuals_meta(1)= seedbank;
        Histo_individuals_meta(2)= to_first_year;
        Histo_individuals_meta(3) = to_second_year;
        Histo_individuals_meta(4) = Histo_individuals_meta(4)+ to_adult_loop;
        
        % Note individuals at May of year +1
        
        %% Disturbance on Adults
        surviving_ratio = (1-sensi_para.proportion_disturbance_Adult);
        Histo_individuals_meta(4:histo_length) = Histo_individuals_meta(4:histo_length).* surviving_ratio;
        
        %% Minimum density correction
        % make sure populations are not infitive small, using all adults as
        % one main class
        Indi_meta_pop = sum(Histo_individuals_meta(4:histo_length));
        Indi_meta_pop_juv = sum(Histo_individuals_meta(2:3));
        if Indi_meta_pop < sensi_para.threshold
            if Indi_meta_pop > adult_test(year,population)*sensi_para.threshold
                histo_tmp = ((Histo_individuals_meta(4:histo_length)./Indi_meta_pop)).*sensi_para.threshold;
                Histo_individuals_meta(4:histo_length) = histo_tmp;
                Histo_stop = 0;
            else
                Histo_individuals(:,population) = Histo_individuals_meta';
                Histo_individuals_meta(4:histo_length) = 0;
                Histo_stop = 1;
            end
            Indi_meta_pop = sum(Histo_individuals_meta(4:histo_length));
        end

        % Make sure there are not infitive small number of seeds and
        % seedlings
        if Histo_individuals_meta(1) < sensi_para.threshold
            if Histo_individuals_meta(1) > seed_test(year,population)*sensi_para.threshold
                Histo_individuals_meta(1) = sensi_para.threshold;
            else
                Histo_individuals_meta(1) = 0;
            end
        end
        if Histo_individuals_meta(2) < sensi_para.threshold
            if Histo_individuals_meta(2) >  class1_test(year,population)*sensi_para.threshold
                Histo_individuals_meta(2) = sensi_para.threshold;
            else
                Histo_individuals_meta(2) = 0;
            end
        end
        if Histo_individuals_meta(3) < sensi_para.threshold
            if Histo_individuals_meta(3) > class2_test(year,population)*sensi_para.threshold
                Histo_individuals_meta(3) = sensi_para.threshold;
            else
                Histo_individuals_meta(3) = 0;
            end
        end
         Indi_meta_pop_juv = sum(Histo_individuals_meta(2:3));
        %% Maximum density correction and lambda calculation
        Individuals_all = sum((Histo_individuals_meta(2:histo_length)));
        % Calcuate lambda and mortality vs establishment
        if  Individuals_all > 0
            % correct for maximum density
            if  Indi_meta_pop >= max_density(1)
                histo_tmp = ((Histo_individuals_meta(4:histo_length)./Indi_meta_pop)).*max_density(1);
                Histo_individuals_meta(4:histo_length) = histo_tmp;
                clear histo_tmp
                Indi_meta_pop = sum(Histo_individuals_meta(4:histo_length));
                if reached == 0                 % !!! add a (population) here but first test the effect of it
                    % put lambda calcuation here, since it will fail the if
                    % loop
                    lambda(population) = (Indi_meta_pop./ (Start_individuals(population))) ^(1/((year+1)-1)) ;
                end
                reached(population) = 1;
            elseif Indi_meta_pop < max_density(2) && reached(population) == 1
                reached(population) = 0; % reset the reached
            end
            
            % calcuate lambda and re-establish populations, note the two
            % year delay to get adult plants (difference between all
            % individuals (classes 2 till end) and Indi_meta_pop (classes 4
            % till end)
            if mortality(population)  == 0
                if   Indi_meta_pop > 0
                    if reached(population) == 0
                        lambda(population) = (Indi_meta_pop./ (Start_individuals(population))) ^(1/((year+1)-1)) ;
                    end
                    % stable state distribution
                end
                Colonisation(population) = 0;
            elseif mortality(population)  == 1
                if  Indi_meta_pop > 0
                    mortality(population) = 0;
                    Colonisation(population) = 1;
                    Start_individuals(population) = sum(Histo_individuals_meta(4:histo_length));
                end
            end
            
        elseif Individuals_all == 0 &&  mortality(population) == 0; % Dissapearing population
            mortality(population) = 1;
            Extinction(population) = 1;
            %else Individuals_all == 0 &&  mortality(population) == 1;
            %Nothing changes
        end
        if Histo_stop == 0
            Histo_individuals(:,population) = Histo_individuals_meta';
        end
        Indi_meta(population) = Indi_meta_pop;
        Indi_meta_juv(population) = Indi_meta_pop_juv;
        Indi_seeds(population) = Histo_individuals_meta(1);
    end % POPULATIONS
    %% Metapopulation growth rates
    total_individual(year+1) = sum(Indi_meta) ;
    Individuals(year+1,:) = Indi_meta;
    Individuals_juv(year+1,:) = Indi_meta_juv; %#ok<NASGU>
    Individuals_seed(year+1,:) = Indi_seeds; %#ok<NASGU>
    Occupancy_adult(year+1) = (length(find(Indi_meta > 0)))./pop_max;
    Occupancy_seeds(year+1) = (length(find(Indi_seeds > 0)))./pop_max;
    if total_individual(year+1) > 0 && (Occupancy_adult(year+1) < 0.95 || pop_max == 3)
        occup_meta =  (Occupancy_adult(year+1)./Occupancy_adult(1))^(1/((year+1)-1));
        %display(lambda_meta) %#ok<NASGU>
        %     else
        %         display('wrong')
    end
    if total_individual(year+1) > 0 && (total_individual(year+1) <max_density(3) || pop_max == 3)
        if (year+1)>= 4
            lambda_meta = (sum(Indi_meta)./sum(Individuals(4,:)))^(1/((year+1)-5));
        else
            lambda_meta = (sum(Indi_meta)./sum(Start_individuals)) ^(1/((year+1)-1)) ;
        end
        if isinf(lambda_meta) ~= 0
            lambda_meta = (sum(Indi_meta)./sum(Start_individuals)) ^(1/((year+1)-1)) ;
        end
    end
    if total_individual(year+1) > 0 && (total_individual(year+1) <max_density(4) || pop_max == 3)
        if (year+1)>= 4
            lambda_meta2 = (sum(Indi_meta)./sum(Individuals(4,:)))^(1/((year+1)-5));
        else
            lambda_meta2 = (sum(Indi_meta)./sum(Start_individuals)) ^(1/((year+1)-1)) ;
        end
        if isinf(lambda_meta) ~= 0
            lambda_meta2 = (sum(Indi_meta)./sum(Start_individuals)) ^(1/((year+1)-1)) ;
        end
    end
    if total_individual(year+1) > 0
        if (year+1)>= 4
            lambda_metaY(year+1) = (sum(Indi_meta)./sum(Individuals(year,:)));
        else
            lambda_metaY(year+1) = NaN;
        end
        if isinf(lambda_metaY(year+1)) ~= 0
            lambda_metaY(year+1) = NaN;
        end
    end
    if total_individual(year+1) > 0
        if (year+1)>= 4
            lambda_metaYcontinue(year+1) =  (sum(Indi_meta)./sum(Individuals(4,:)))^(1/((year+1)-5));
        else
            lambda_metaYcontinue(year+1) = NaN;
        end
        if isinf(lambda_metaYcontinue(year+1)) ~= 0
            lambda_metaY(year+1) = NaN;
        end
    end
    
    clear Seeds_total
    %% Extinction and Colonisation rates
    if year == 1
        Empty_rate(year) = sensi_para.zero_densities./pop_max;
        Colonisation_rate(year) = (1-sensi_para.zero_densities);
    end
    Empty_rate(year+1) = (sum(mortality))./pop_max;
    Colonisation_rate(year+1) = (sum(Colonisation))./pop_max;
    Extinction_rate(year+1) = (sum(Extinction))./pop_max; %#ok<NASGU>
    
end % YEARS
cd temp
name_file = ['output','_',int2str(run),'.mat'];
save(name_file,'lambda','Individuals','Individuals_juv',...
    'Histo_individuals','Occupancy_adult','Occupancy_seeds','Individuals_seed',...
    'Start_individuals','lambda_meta','Existing_Stable_state','Prop_among',...
    'sensi_para','total_individual', 'occup_meta','Empty_rate', 'Colonisation_rate',...
    'Extinction_rate','lambda_meta2','lambda_metaY','lambda_metaYcontinue')
cd ..
% end run