% Stochastic IPM Base model for a connected B. olerecea metapopulation,
% including dispersal

% Outputs: lambda and stable state distribution per population

% The model simulates population size per year in a point population with a
% diameter of 50-meters. In between populations is not habitable habitat

% the model is based on number of plants flowering or not before seed set.
%The model starts with seed set and incorperation of the % plants flowering per class
tic
clc
clear all
warning('off')
mkdir('temp')
clear all
clear all

% General parameters
Sensis = [1:19];
definition_function = 'name.m';
extention = 'Extention';

folder = 'definitions/';
definition_function = [folder,definition_function];
folder_out = 'Metpop_outcome_runs/';
if exist(folder_out) == 0
    mkdir(folder_out)
end
if exist('Output.mat') ~= 0
    copyfile ('Output.mat', 'OutputBK.mat')
    load('Output.mat')
    %delete('Output.mat')
else
    Output.org = {[1,2]};
end

if exist('define_function') ~= 0
    delete 'define_function.m'
end
copyfile (definition_function, 'define_function.m')
clear definition_function

max_diameter = 100;
histo_length =  max_diameter +3;
size_steps = ceil(max_diameter/10);

cd Data_analyses_results
load('pops&dispersal.mat')
cd ../
max_density(1) = (max(max(Dorset_populations.densities))).*2; % maximum number of adult plants
max_density(2) = 0.5.*max_density(1); % max density tolerance

%save('parameters_start','Sensis','max_diameter','histo_length','max_density');
% Define all populations IPM paeameters
Oleracea_parameters_load;
%load('parameters_start.mat')
population_types = Oleracea_population_types;
save('start_parameters')


Time_run = 0;
pops_only = 0;
%% Run the model for every sensi provided
for sensis = 1:1:length(Sensis)
    sensi = Sensis(sensis);
    tic;
    save('Output','Output')
    run = 1;
    define_function
    max_density(3) = (max_density(1).*0.5).*pop_max;   
    max_density(4) = (max_density(1).*0.95).*pop_max; 
    delete('Prop_among.mat')
    if pop_max == 3
        pops_only = 1;
    end
    
    %% Calculate Kernels
    Oleracea_Seed_flow_kernels
    toc
    save('Prop_among', 'Prop_among');
    
    %% Send runs to processor
    %     number_of_processors = getenv('NUMBER_OF_PROCESSORS');
    %     number_of_processors = str2double(number_of_processors);
    stints = floor(run_max/(10-1));
    output = 1;
    if stints > 0
        display('Paralel job')
        job = createJob('configuration', 'Full');
        for run = 1:run_max
            %display(run)
            createTask(job, @Oleracea_actual_simulations, 1,...
                {run,population_types,max_density,Prop_among,Sensis,...
                max_diameter,size_steps,sensi,histo_length,run_max,extention});
        end
        display(' ')
        display('Running Sensi:')
        disp(sensi)
        submit(job);
        waitForState(job, 'finished');
        results = getAllOutputArguments(job);
        destroy(job)
    else
        display('No paralel')
        for run = 1:run_max
            display(' ')
            display('Running Sensi:')
            disp(sensi)
            [run_max] = Oleracea_actual_simulations...
                (run,population_types,max_density,Prop_among,Sensis,...
                max_diameter,size_steps,sensi,histo_length,run_max,extention);
        end
    end
    
    %% Collect output per run and do statistics and write to Output file
    load('Output.mat')
    Output.population_types = population_types;
    if sensi ~= sensi_para.baseline
      Lambda_meta_base = Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_meta_base;
      Occup_meta_base=  Output.(genvarname(['s',int2str(sensi_para.baseline)])).Occup_meta_base;
      Lambda_base =  Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_base;
    end
    
    
    cd temp
    for run = 1:run_max
        name_file = ['output','_',int2str(run),'.mat'];
        load(name_file)
        Lambda(run, :) =  lambda; %#ok<*SAGROW>
        Individualss(run,:,:) = Individuals;
        Individualss_seeds(run, :,:) = Individuals_seed;
        Individualss_juv(run, :,:) = Individuals_juv;
        Histo_individuals_metas(run,:,:) = Histo_individuals;
        Occupancies(run,:) = Occupancy_adult;
        Occupancies_seeds(run,:) = Occupancy_seeds;
        Start_individualss(run,:) = Start_individuals;
        Lambda_metaY(run,:) = lambda_metaY;
        Lambda_meta(run) = lambda_meta;
        Lambda_meta2(run) = lambda_meta2;
        Occup_meta(run) = occup_meta;
        Existing_Stable_states(run,:,:) = Existing_Stable_state;
        Total_individuals(run,:) = total_individual;
        Empty_rates(run,:) = Empty_rate;
        Colonisation_rates(run,:) = Colonisation_rate;
        Extinction_rates(run,:) = Extinction_rate;
        
        if exist('Do_sensitivity') ~= 0
            if strcmpi('y',Do_sensitivity) ~= 0
                % sensitivity per run for meta population lambda
                if sensi ~= sensi_para.baseline
                    base = Lambda_meta_base(run);% Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_meta.mean;
                    new = lambda_meta;
                    if isnan(base) ~= 1 && isnan(new) ~= 1 % This is SENSITIVITY !!!
                        %Lambda_meta
                        sens_tmp = (new/base) -1;
                        if new < base
                            sens_tmp = 1-(1./(base/new));
                        end
                        if sensi_para.change_factor ~=1
                            sens_meta(run) = sens_tmp./abs(1-sensi_para.change_factor);
                        else
                            sens_meta(run) = 0;
                        end
                    end
                else
                    sens_meta(run) = 0;
                end
                clear sens_tmp
                
                %% sensitivity per run for meta occupancy
                if sensi ~= sensi_para.baseline
                    base = Occup_meta_base(run);% Output.(genvarname(['s',int2str(sensi_para.baseline)])).Occupancy_meta.occupancy(1);
                    new = occup_meta;
                    if isnan(base) ~= 1 && isnan(new) ~= 1 % This is SENSITIVITY !!!
                        %Lambda_meta
                        sens_tmp = (new/base) -1;
                        if new < base
                            sens_tmp = 1-(1./(base/new));
                        end
                        if sensi_para.change_factor ~=1
                            occu_meta(run) = sens_tmp./abs(1-sensi_para.change_factor);
                        else
                            occu_meta(run) = 0;
                        end
                    end
                else
                    occu_meta(run) = 0;
                end
                clear sens_tmp
                
                
                %% sensitivity per run for populations lambda
                value = lambda;
                value(isnan(value)==1) = [];
                new =  mean(mean(value));
                clear value
                if sensi ~= sensi_para.baseline
                    value = Lambda_base(run,:); %Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_mean.mean;
                    %Output.Lambda(1,:,sensi_para.baseline);
                    value(isnan(value)==1) = [];
                    base =  mean(mean(value));
                    if isnan(base) ~= 1 && isnan(new) ~= 1 % This is SENSITIVITY !!!
                        %Lambda_meta
                        sens_tmp = (new/base) -1;
                        if new < base
                            sens_tmp = 1-(1./(base/new));
                        end
                        if sensi_para.change_factor ~=1
                            sens_pop(run) =  sens_tmp./abs(1-sensi_para.change_factor);
                        else
                            sens_pop(run) = 0;
                        end
                    end
                else
                    sens_pop(run) = 0;
                end
                clear sens_tmp

                if strcmpi('y',Do_sensitivity) == 0
                    sensi_para.change_factor = 1;
                end
            end
        end
        
        diror = ['../','Store_',extention,'/','outputs','_',int2str(sensi)];
        if exist(diror) == 0
            mkdir(diror)
        end
        
        dir2 =  [diror,'/'];
        %name_files = ['temp/','sensitivity_factors','_',int2str(run),'.mat'];
        movefile(name_file,dir2);
        clear dir2
        clear diror
        clear lambda
        clear Individuals
        clear Individuals_seed
        clear Histo_individual
        clear Occupancy_adult
        clear Occupancy_seeds
        clear Start_individuals
        clear lambda_meta;
        clear lambda_meta2;
        clear lambda_metaY;
        clear ocup_meta;
        clear Existing_Stable_state
        clear sens_tmp
        clear value
        clear new
        clear new_std
        clear base
        clear occup_meta
        clear total_individual
        clear Empty_rate
        clear Colonisation_rate
        clear Extinction_rate
    end % per run
    % save baseline
    
    cd ..
    if sensi == sensi_para.baseline
        Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_meta_base = Lambda_meta;
        Output.(genvarname(['s',int2str(sensi_para.baseline)])).Occup_meta_base = Occup_meta;
        Output.(genvarname(['s',int2str(sensi_para.baseline)])).Lambda_base = Lambda;
    end
    
    %% Statistics and Save output
    
    
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.mean = (mean(Lambda_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.std = (std(Lambda_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.median = (median(Lambda_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.per5 = (prctile(Lambda_meta,5))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.per25 = (prctile(Lambda_meta,25))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.per75 = (prctile(Lambda_meta,75))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.per95 = (prctile(Lambda_meta,95))';
    
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.mean = (mean(Lambda_meta2))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.std = (std(Lambda_meta2))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.median = (median(Lambda_meta2))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.per5 = (prctile(Lambda_meta2,5))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.per25 = (prctile(Lambda_meta2,25))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.per75 = (prctile(Lambda_meta2,75))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta95.per95 = (prctile(Lambda_meta2,95))';
    
    if exist('Do_sensitivity') ~= 0
        if strcmpi('y',Do_sensitivity) ~= 0
            Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.sensitivity(1) = mean(sens_meta);
            Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.sensitivity(2) = std(sens_meta);
            Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.sensitivity(1) = mean(occu_meta);
            Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.sensitivity(2) = std(occu_meta);
        end
    end
    
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(1) = (mean(Occup_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(2) = (std(Occup_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(3) = (median(Occup_meta))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(4) = (prctile(Occup_meta,5))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(5) = (prctile(Occup_meta,25))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(6) = (prctile(Occup_meta,75))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.occupancy(7) = (prctile(Occup_meta,95))';
    
     if exist('Do_sensitivity') ~= 0
        if strcmpi('y',Do_sensitivity) ~= 0
                Output.(genvarname(['s',int2str(sensi)])).Lambda_pop.sensitivity(1)  = mean(sens_pop);
        Output.(genvarname(['s',int2str(sensi)])).Lambda_pop.sensitivity(2) = std(sens_pop);
        end
     end
     
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.mean = (mean(Lambda_metaY(:,5:year_max),1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.std = (std(Lambda_metaY(:,5:year_max),0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.median = (median(Lambda_metaY(:,5:year_max),1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.per5 = (prctile(Lambda_metaY(:,5:year_max),5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.per25 = (prctile(Lambda_metaY(:,5:year_max),25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.per75 = (prctile(Lambda_metaY(:,5:year_max),75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.per95 = (prctile(Lambda_metaY(:,5:year_max),95,1))';
    Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.meanPeriod = mean(Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.mean);
     Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.MedianPeriod = median(Output.(genvarname(['s',int2str(sensi)])).Lambda_meta_Year.mean);
    
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.mean = (mean(Total_individuals,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.std = (std(Total_individuals,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.median = (median(Total_individuals,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.per5 = (prctile(Total_individuals,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.per25 = (prctile(Total_individuals,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.per75 = (prctile(Total_individuals,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Total_Individuals.per95 = (prctile(Total_individuals,95,1))';
    
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.mean = (mean(Occupancies,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.std = (std(Occupancies,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.median = (median(Occupancies,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.per5 = (prctile(Occupancies,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.per25 = (prctile(Occupancies,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.per75 = (prctile(Occupancies,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.adult.per95 = (prctile(Occupancies,95,1))';
    
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.mean = (mean(Occupancies_seeds,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.std = (std(Occupancies_seeds,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.median = (median(Occupancies_seeds,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.per5 = (prctile(Occupancies_seeds,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.per25 = (prctile(Occupancies_seeds,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.per75 = (prctile(Occupancies_seeds,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Occupancy.seeds.per95 = (prctile(Occupancies_seeds,95,1))';
    
    Output.(genvarname(['s',int2str(sensi)])).Individuals.mean = squeeze(mean(Individualss,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.std = squeeze(std(Individualss,0,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.median = squeeze(median(Individualss,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.per5 = squeeze(prctile(Individualss,5,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.per25 = squeeze(prctile(Individualss,25,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.per75 = squeeze(prctile(Individualss,75,1));
    Output.(genvarname(['s',int2str(sensi)])).Individuals.per95 = squeeze(prctile(Individualss,95,1));
    
    Output.(genvarname(['s',int2str(sensi)])).Seeds.mean = squeeze(mean(Individualss_seeds,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.std = squeeze(std(Individualss_seeds,0,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.median = squeeze(median(Individualss_seeds,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.per5 = squeeze(prctile(Individualss_seeds,5,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.per25 = squeeze(prctile(Individualss_seeds,25,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.per75 = squeeze(prctile(Individualss_seeds,75,1));
    Output.(genvarname(['s',int2str(sensi)])).Seeds.per95 = squeeze(prctile(Individualss_seeds,95,1));
    
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.mean = squeeze(mean(Individualss_juv,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.std = squeeze(std(Individualss_juv,0,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.median = squeeze(median(Individualss_juv,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.per5 = squeeze(prctile(Individualss_juv,5,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.per25 = squeeze(prctile(Individualss_juv,25,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.per75 = squeeze(prctile(Individualss_juv,75,1));
    Output.(genvarname(['s',int2str(sensi)])).Juveniles.per95 = squeeze(prctile(Individualss_juv,95,1));
    
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.mean = (mean(Colonisation_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.std = (std(Colonisation_rates,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.median = (median(Colonisation_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.per5 = (prctile(Colonisation_rates,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.per25 = (prctile(Colonisation_rates,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.per75 = (prctile(Colonisation_rates,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Colonisation.per95 = (prctile(Colonisation_rates,95,1))';
    
    Output.(genvarname(['s',int2str(sensi)])).Extinction.mean = (mean(Extinction_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.std = (std(Extinction_rates,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.median = (median(Extinction_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.per5 = (prctile(Extinction_rates,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.per25 = (prctile(Extinction_rates,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.per75 = (prctile(Extinction_rates,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Extinction.per95 = (prctile(Extinction_rates,95,1))';
    
    Output.(genvarname(['s',int2str(sensi)])).Empty.mean = (mean(Empty_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.std = (std(Empty_rates,0,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.median = (median(Empty_rates,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.per5 = (prctile(Empty_rates,5,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.per25 = (prctile(Empty_rates,25,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.per75 = (prctile(Empty_rates,75,1))';
    Output.(genvarname(['s',int2str(sensi)])).Empty.per95 = (prctile(Empty_rates,95,1))';
    
    
    %% Stable states
    for i = 2:2:(histo_length-3)
        Stable_state_runs(:,(i/2),:) = (sum(Histo_individuals_metas(:,(i+2):(i+3),:),2))./(sum(Histo_individuals_metas(:,4:histo_length,:),2));
    end
    Stable_state_S_ratio(:,1,:) = (sum(Histo_individuals_metas(:,2,:),2))./(sum(Histo_individuals_metas(:,4:histo_length,:),2));
    Stable_state_S_ratio(:,2,:) = (sum(Histo_individuals_metas(:,3,:),2))./(sum(Histo_individuals_metas(:,4:histo_length,:),2));
    if pop_max <= 10
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.existing.Mean = squeeze(mean(Existing_Stable_states));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.mean = squeeze(mean(Stable_state_runs,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.median = squeeze(median(Stable_state_runs,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.std = squeeze(std(Stable_state_runs,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.per5 = squeeze(prctile(Stable_state_runs,5,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.per25 = squeeze(prctile(Stable_state_runs,25,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.per75 = squeeze(prctile(Stable_state_runs,75,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.per95 = squeeze(prctile(Stable_state_runs,95,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.Sratio.Mean = squeeze(mean(Stable_state_S_ratio,1));
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.Sratio.std = squeeze(std(Stable_state_S_ratio,1));
    else
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.existing.Mean = mean(squeeze(mean(Existing_Stable_states)),2);
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.mean = mean(squeeze(mean(Stable_state_runs,1)),2);
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.median = median(squeeze(median(Stable_state_runs,1)),2);
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.std = mean(squeeze(std(Stable_state_runs,1)),2);
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.Sratio.Mean = mean(squeeze(mean(Stable_state_S_ratio,1)),2);
        Output.(genvarname(['s',int2str(sensi)])).Stable_state.Sratio.std = mean(squeeze(std(Stable_state_S_ratio,1)),2);
        
        %         if run_max <=10
        %             Output.(genvarname(['s',int2str(sensi)])).runs.Lambda = Lambda;
        %             Output.(genvarname(['s',int2str(sensi)])).runs.Individuals = Individualss;
        %             Output.(genvarname(['s',int2str(sensi)])).runs.Stable_state = Stable_state_runs;
        %             Output.(genvarname(['s',int2str(sensi)])).runs.Start_individuals = Start_individualss;
        %             Output.(genvarname(['s',int2str(sensi)])).runs.lambda_meta = Lambda_meta;
        %         end
    end
    clear Individualss
    clear Stable_state_runs
    clear Individualss_seeds
    clear Individualss_juv
    clear Histo_individuals_metas
    clear Occupancies
    clear Occupancies_seeds
    clear Start_individualss
    clear Lambda_meta
    clear Existing_Stable_states
    clear Total_individuals
    clear Stable_state_runs
    clear Colonisation_rates
    clear Extinction_rates
    clear Empty_rates
    clear Lambda_metaY
    clear Lambda_meta2
    clear Occup_meta
    clear Start_individualss
    clear  Existing_Stable_states
    clear Stable_state_S_ratio
    clear Lambda_meta_base
    clear Occup_meta_base
    clear Lambda_base
    clear results
    clear prop_disp
    %% Population Lambdas
    for population = 1:1:pop_max
        Lambda_pop = Lambda(:,population);
        Lambda_pop(Lambda_pop ==0) = [];
        if isempty(Lambda_pop) ~= 1
            Lambdas(1,population) = mean(Lambda_pop);
            Lambdas(2,population) = median(Lambda_pop);
            Lambdas(3,population) = std(Lambda_pop);
            Lambdas(4,population) =  prctile(Lambda_pop,5);
            Lambdas(5,population) =  prctile(Lambda_pop,25);
            Lambdas(6,population) =  prctile(Lambda_pop,75);
            Lambdas(7,population) =  prctile(Lambda_pop,95);
        else
            Lambdas(1:7,population) = NaN;
        end
        clear Lambda_pop
    end   
    clear Lambda
    
    value =  Lambdas(1,:);
    value(isnan(value)==1) = [];
    Output.(genvarname(['s',int2str(sensi)])).Lambda_mean.mean =  mean(value);
    Output.(genvarname(['s',int2str(sensi)])).Lambda_mean.std = std(value);
    if pop_max <= 100
        Output.Lambda(:,:,sensi) = Lambdas;
    else
        clear Lambdas
    end
    
    %% Store parameters used
    
    Output.(genvarname(['s',int2str(sensi)])).kernel = kernel;
    Output.(genvarname(['s',int2str(sensi)])).Prop_among = Prop_among;
    delete('Prop_among.mat')
    
    sensi_para = rmfield(sensi_para,'show');
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para = sensi_para;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.pop_max = pop_max;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.run_max = run_max;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.year_max = year_max;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.max_diameter = max_diameter;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.size_steps = size_steps;
    Output.(genvarname(['s',int2str(sensi)])).Sensi_para.max_density = max_density;
    Output.Sensi_para_baseline = sensi_para.baseline;
    
    if exist('Do_sensitivity') ~= 0
        Output.(genvarname(['s',int2str(sensi)])).Sensi_para.Do_sensitivity = Do_sensitivity;
    end
    if exist('zero_densities')~= 0
        Output.(genvarname(['s',int2str(sensi)])).Sensi_para.zero_densities = zero_densities;
    end
    
    %% Screen output per sensi
    if  pops_only == 1
        to_display = Output.(genvarname(['s',int2str(sensi)])).Lambda_mean.mean;
        display('Average per population lambda')
        disp(to_display)
        clear to_display
        
            
     if exist('Do_sensitivity') ~= 0
        if strcmpi('y',Do_sensitivity) ~= 0
        to_display = Output.(genvarname(['s',int2str(sensi)])).Lambda_pop.sensitivity(1);
        display('Sensitivity:')
        disp(to_display)
        end
     end
    else
        to_display = Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.mean;
        display('Meta population lambda')
        disp(to_display)
        
        if exist('Do_sensitivity') ~= 0
            if strcmpi('y',Do_sensitivity) ~= 0
                to_display = Output.(genvarname(['s',int2str(sensi)])).Lambda_meta50.sensitivity(1);
                display('Sensitivity:')
                disp(to_display)
                
                to_display = Output.(genvarname(['s',int2str(sensi)])).Occupancy_meta.sensitivity(1);
                display('Occupancy Sensitivity:')
                disp(to_display)
            end
        end
    end
    
    
    sensi_toc = toc;
    Output.(genvarname(['s',int2str(sensi)])).time_run = toc;
    toc
    Time_run  = Time_run + toc;
    Output.Time_run = Time_run;
    clear Mortality
    clear Mortality_year
    clear Established
    clear Lambda
    clear kernels
    clear Prop_amongs
    clear sensi_para
end % SENSIS

%% Create Screen output after finishing
Series = [ceil(year_max/20),ceil(year_max/10),ceil(year_max/4),ceil(year_max/2),year_max];
for i = 1:length(Sensis)
    Output.Lambdas_meta50(i,2) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.mean;
    Output.Lambdas_meta50(i,3) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.std;
    Output.Lambdas_meta50(i,4) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.per5;
    Output.Lambdas_meta50(i,5) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.per95;
    
    Output.Lambdas_meta95(i,2) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.mean;
    Output.Lambdas_meta95(i,3) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.std;
    Output.Lambdas_meta95(i,4) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.per5;
    Output.Lambdas_meta95(i,5) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.per95;
    
    Output.Lambdas_mean(i,2) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_mean.mean;
    Output.Lambdas_mean(i,3) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_mean.std;

    Output.Occupancy_meta(i,2) = Output.(genvarname(['s',int2str(Sensis(i))])).Occupancy_meta.occupancy(1);
    Output.Occupancy_meta(i,3) = Output.(genvarname(['s',int2str(Sensis(i))])).Occupancy_meta.occupancy(2);

    
    if exist('Do_sensitivity') ~= 0
        if strcmpi('y',Do_sensitivity) ~= 0
            Output.Lambdas_meta50(i,6) =Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.sensitivity(1);
            Output.Lambdas_meta50(i,7) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_meta50.sensitivity(2);
            Output.Lambdas_mean(i,4) =Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_pop.sensitivity(1);
            Output.Lambdas_mean(i,5) = Output.(genvarname(['s',int2str(Sensis(i))])).Lambda_pop.sensitivity(2);
            Output.Occupancy_meta(i,4) = Output.(genvarname(['s',int2str(Sensis(i))])).Occupancy_meta.sensitivity(1);
            Output.Occupancy_meta(i,5) = Output.(genvarname(['s',int2str(Sensis(i))])).Occupancy_meta.sensitivity(2);
        end
    end

    if pops_only == 1
        Output.Lambdas.mean(i,:) =  Output.Lambda(1,:,i);
        Output.Lambdas.std(i,:) =  Output.Lambda(3,:,i);
    end
end
Output.Lambdas_mean(:,1) = Sensis';
Output.Lambdas_meta50(:,1) = Sensis';
Output.Lambdas_meta95(:,1) = Sensis';
Output.Occupancy_meta(:,1) = Sensis';

if pops_only == 1
   for i = 1:length(Sensis)
    Output.Lambdas_pops.mean(i,1) = Output.Lambda(1,1,i);
    Output.Lambdas_pops.mean(i,2) = Output.Lambda(1,2,i);
    Output.Lambdas_pops.mean(i,3) = Output.Lambda(1,3,i);
       
    Output.Lambdas_pops.std(i,1) = Output.Lambda(3,1,i);
    Output.Lambdas_pops.std(i,2) = Output.Lambda(3,2,i);
    Output.Lambdas_pops.std(i,3) = Output.Lambda(3,3,i);
   end
end
    


% sensitivity of occurences
%Oleracea_occupancy_sens;
save('Output','Output')
delete('parameters.mat')
delete ('parameters_start.mat')
delete('start_parameters.mat')
rmdir('temp','s')

name_file = [folder_out,'Output','_',extention,'.mat'];
movefile ('Output.mat', name_file)
save('extention','extention','folder_out', 'pops_only')

%% Clear and show screen output
delete('define_function.m')
delete('sensi_tmp.mat')

clear all
clc
clear all
clc
clear all

load('extention.mat')
delete('extention.mat')
name_file = [folder_out,'Output','_',extention,'.mat'];
load(name_file)

if pops_only == 1
    display(Output.Lambdas_pops.mean)
else
    display(Output.Lambdas_meta50)
    display(Output.Occupancy_meta)
end



Seconds_run = Output.Time_run
Output = rmfield(Output,'org');
clear Seconds_run
clear extention
clear name_file
clear folder_out
clear pops_only