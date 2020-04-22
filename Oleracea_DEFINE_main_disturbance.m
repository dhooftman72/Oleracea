% run and pop parameters

%% Set Base variables
pop_max = 1000;
run_max = 1;
year_max = 100;
sensi_para.area_size = 12;
sensi_para.area_size_sens_wind = sensi_para.area_size;
sensi_para.area_size_sens_human = sensi_para.area_size;
sensi_para_integral_steps = 1;
sensi_para.wind_dispersal_present = 'Y';
sensi_para.human_dispersal_present = 'Y';
sensi_para.pick_up_sensitivity = 1;
sensi_para.threshold = 1;
sensi_para.Pu = 0.023;
sensi_para.show = ceil(100./pop_max);
sensi_para.proportion_disturbance = 0.25; % add here the proportion disturbance per site
sensi_para.all_years = 'N'; % all years together (yes) or draw from single years (N)
sensi_para.proportion_disturbance_Adult = 0;
sensi_para.mortality_disturbance_Adult = 0;
sensi_para.zero_densities = 1-(0.35);
sensi_para.baseline = 1; % baseline run for sensitivity calculations
Do_sensitivity = 'N';
sensi_para.change_factor = 1;

% % set histogram diagram
Histo_base = zeros(max_diameter,pop_max);
Histo_sizes(:,1) = 2:(max_diameter+1);

%% Set distances

% create 1D pop_max populations, at area_size distance apart, starting with 0.
Distances_centre = 0:sensi_para.area_size:(sensi_para.area_size*(pop_max-1));
% Distances_edge(1) = 0;
% Distances_edge_2 = (sensi_para.area_size/2):sensi_para.area_size:(sensi_para.area_size*(pop_max-1));
% Distances_edge(2:(length(Distances_edge_2)+1)) = Distances_edge_2 ;
clear Distances_edge_2
for i = 1:pop_max
    for y = 1:1:pop_max
        Distance_matrix.centre(i,y)= abs(Distances_centre(i)-Distances_centre(y));
%         Distance_matrix.edge(i,y)= abs(Distances_edge(i)-Distances_edge(y));
    end
end
Distance_matrix.oneD.centre = Distances_centre;
% Distance_matrix.oneD.edge = Distances_edge;
clear Distances_edge
clear Distances_centre
clear i
clear y



%% create populations
cd Results
load('pops&dispersal.mat')
cd ../
mean_densities = mean(Dorset_populations.densities,2);
mean_densities = mean_densities./(prctile(mean_densities,95));
mean_densities(mean_densities>1) = 1;
histo = hist(mean_densities(mean_densities>0),size_steps);
histo(2,:) = 1/size_steps:(1/size_steps):1;
pdf_pops = histo(1,:)/sum(histo(1,:));
cdf_pops = cumsum(pdf_pops);

if sensi ~= sensi_para.baseline
    cd temp
    name_file = ['sensitivity_factors','_',int2str(run),'.mat'];
    load(name_file)
    cd ..
end

count_zero = 0;
again = 1;
while again ~= 0
    for population = 1:1:pop_max
        
        type_draw = ceil(rand*3);
        if sensi ~= sensi_para.baseline
            type_draw = types_draw(population);
        else
            types_draw(population) = type_draw; %#ok<*SAGROW>
        end
         hist_draw =  min(find(cdf_pops>= rand)); %#ok<*MXFND>
        if sensi ~= sensi_para.baseline
            hist_draw = hist_draws(population);
        else
            hist_draws(population) =  hist_draw;
        end
        
        zero_draw = rand;
        if sensi ~= sensi_para.baseline
            zero_draw = zero_draws(population);
        else
            zero_draws(population) =  zero_draw;
        end
        
        if type_draw == 1
            populations.(genvarname(['population','_',int2str(population)])) = population_types.Kim;
        elseif type_draw == 2
            populations.(genvarname(['population','_',int2str(population)])) = population_types.OH;
        elseif type_draw == 3
            populations.(genvarname(['population','_',int2str(population)])) = population_types.Win;
        else
            display('error, unknown population type')
        end
        if zero_draw <= sensi_para.zero_densities
            populations.(genvarname(['population','_',int2str(population)])).Histo_base =...
                populations.(genvarname(['population','_',int2str(population)])).Histo_base.*0;
            count_zero = count_zero +1;
        else
            populations.(genvarname(['population','_',int2str(population)])).Histo_base = ...
                populations.(genvarname(['population','_',int2str(population)])).Histo_base.*histo(2,hist_draw);
        end
    end
    zero_rate =  count_zero./pop_max;
    
    if  zero_rate  < (sensi_para.zero_densities - (1/pop_max)) ||  zero_rate > (sensi_para.zero_densities + (1/pop_max)) || zero_rate == 1
        again = 1;
        zero_rate = 0;
        count_zero = 0;
    else
        again = 0;
    end
    if zero_rate == 1 && again == 0
        display('Something went terribly wrong')
        ccc
    end
end

%% set changed parameters
if sensi == 1
    sensi_para.proportion_disturbance = 1; % add here the proportion disturbance per site
    display('Disturbed')
elseif sensi == 2
    sensi_para.proportion_disturbance = 0.9;
elseif sensi == 3
    sensi_para.proportion_disturbance = 0.85;    
elseif sensi == 4
    sensi_para.proportion_disturbance = 0.8;
elseif sensi == 5
    sensi_para.proportion_disturbance = 0.75;
elseif sensi == 6
    sensi_para.proportion_disturbance = 0.7;
elseif sensi == 7
    sensi_para.proportion_disturbance = 0.65;
elseif sensi == 8
    sensi_para.proportion_disturbance = 0.6;
elseif sensi == 9
    sensi_para.proportion_disturbance = 0.55;    
elseif sensi == 10
    sensi_para.proportion_disturbance = 0.5;
 elseif sensi == 11
    sensi_para.proportion_disturbance = 0.45;   
elseif sensi == 12
    sensi_para.proportion_disturbance = 0.4;
elseif sensi == 13
    sensi_para.proportion_disturbance = 0.35;
elseif sensi == 14
    sensi_para.proportion_disturbance = 0.3;
elseif sensi == 15
    sensi_para.proportion_disturbance = 0.25;
elseif sensi == 16
    sensi_para.proportion_disturbance = 0.2;
elseif sensi == 17
    sensi_para.proportion_disturbance = 0.15;
elseif sensi == 18
    sensi_para.proportion_disturbance = 0.1;
elseif sensi == 19
    sensi_para.proportion_disturbance = 0;
    display('UnDisturbed')
end

if sensi == sensi_para.baseline
    if strcmpi('y',sensi_para.all_years) ~= 0
        population_years (year_max) = 7;
    else
        population_years = randi([1,6],year_max);
    end
    cd temp
    name_file = ['sensitivity_factors','_',int2str(run),'.mat'];
    save(name_file,'population_years', 'types_draw', 'hist_draws','zero_draws')
    cd ..
end