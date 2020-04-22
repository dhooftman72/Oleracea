% in the module, the kernels for seed dispersal by wind and human are
% calcuated,
cd Data_analyses_results
load('pops&dispersal.mat')
cd ../

%% Create the combined wind dispersal kernel for area size intervals
%initate matrix
clear Prop_among
kernel = 0;
max_distance = sensi_para.area_size*(pop_max-1);
%Make unique distance matrices
%Steps are distance centre
step = sensi_para_integral_steps;
size = sensi_para.area_size;
% creating bins from Wichman et al. 2009, following Soons
% et al. 2004
tot_bins = ceil(max(Wind_dispersal_Soons.extreme))./step; % so this is the lower edge;
prop_disp(1:tot_bins) = hist(Wind_dispersal_Soons.extreme,[tot_bins])./50000; 
prop_disp((tot_bins+1):step:(max_distance.*2)) = 0;
kernel.wind(1,:) = prop_disp;
total = length(kernel.wind(1,:));
kernel.wind(2,1:total) = 0:step:(step*(total-1));
%% Human mediated dispersal following Wichman et al. 2009
% create Kernel for size_area intervals
% Make the attacchment into a PDF, with the CDF adding to "a". By
% calcuating which proportion of the orginal drops of during the patch.
% A small proportion beyond the populations dissapears,
% if instead of x + area_size  "inf" is added, the cum_sum of this
% PDF is 'a'.
Human_disp_a = 0.947; % power exponential from Wichman et al 2009
Human_disp_b = 0.165; % power exponential from Wichman et al 2009
fun_num = @(x) Human_disp_a*(exp(-(x.^Human_disp_b))); % the proportion of seeds attached at start patch at distance x
for i = 1:1:(max_distance.*2)
    x = i-step;
    kernel.attached(1,i) = fun_num(x);
    kernel.attached(2,i) = x;
    % proportion still attached at start of patch
end

% Human dispersal is from the edge of the patch to a next patch, where
% it is droppend off at every interval. Hence the drop rate per
% interval is summed as integral with stepsize 'step', but no integral
% over x1 is taken, since only the end is important = size -step
% Within the patch follows the same logic, only seeds are picked off
% and dropped at every intreval in the patch. This is seen as a
% directional movement over x-intervals, resulting in a value still attached at the end
% of the patch.
for i = 1:1:length(Distance_matrix.oneD.centre)
        count = 0;
        x2 = (Distance_matrix.oneD.centre(i)); % beginning of the patch
        x3 = (Distance_matrix.oneD.centre(i)+size)-step; % end of the patch        
        for j = step:step:(size)
            x1 = (Distance_matrix.oneD.centre(1)+size)-j; % distance from source
            dist_1 = find(kernel.attached(2,:)==(x2-x1));
            dist_2 = find(kernel.attached(2,:)==(x3-x1));
            if dist_1 ~= dist_2
                count = count +1;
                disp(count) = (kernel.attached(1,dist_1) - kernel.attached(1,dist_2));
                % Note no correction for the annulus, since this is one
                % directional already
            end
        end
        kernel.dispersal_human(i) = sum(disp)/count;
        clear j
    clear disp
end

%% Wind dispersal following Soons et al. 2004
% as 4-way intergral over every x and every y of source and target.
maximum_distance = ((ceil(max(Wind_dispersal_Soons.extreme)./size)+1));
if length(Distance_matrix.oneD.centre) < maximum_distance 
    maximum_distance = length(Distance_matrix.oneD.centre);
end
for i = 1:1:maximum_distance
    count = 0;
    for x1 = 0:step:((size)-step)
        for y1 = 0:step:((size)-step)
            for x2 = Distance_matrix.oneD.centre(i):step:((Distance_matrix.oneD.centre(i)+size)-step)
                for y2 = 0:step:((size)-step)
                    count = count + 1;
                    distance = sqrt(((x1-x2).^2) + ((y1-y2).^2));
                    distance = (ceil(distance./step)).*step;
                    dist = find(kernel.wind(2,:)==distance);
                    disp(count) = kernel.wind(1,dist)./((size.*(1/step)).^2);    
                    % correction for number of x, y points of source
                    % note not of target, since those need to be summed
                    % over the 144 values there
                    if distance ~= 0
                        disp_annu(count) = disp(count) * (step./((2*pi).*distance)); % correction for annulus
                    end
                end
            end
        end
    end
    clear x1
    clear x2
    clear y1
    clear y2
    kernel.dispersal_wind_annulus(i) =  sum(disp_annu);
    kernel.dispersal_wind_lost(i) =  (sum(disp)./(size*(1/step)))-(2*kernel.dispersal_wind_annulus(i));
    clear disp
    clear disp_annu
end
kernel.dispersal_wind_annulus(i+1:pop_max+1) = 0;
kernel.dispersal_wind_lost(i+1:pop_max+1) = 0;
%% calcute the wind dispersal matrix throughout the meta_population from Kernel
for i = 1:pop_max
    for y = 1:1:pop_max
        x = Distance_matrix.centre(i,y);
        if x > 0
            Prop_among.wind_annulus(i,y) = kernel.dispersal_wind_annulus((find(Distance_matrix.oneD.centre == x)));
            if x == size
                Prop_among.wind_annulus(i,y) = Prop_among.wind_annulus(i,y);
            end
        end         
    end
    Prop_among.wind_lost = sum(kernel.dispersal_wind_lost((1+1):length(kernel.dispersal_wind_lost)));
    clear  wind_edge_1
    clear  wind_edge_2
end
clear x
clear i
clear y

% calculate two sided dissapearences
disapearance_tmp = sum( Prop_among.wind_annulus,1);
Prop_among.wind_edge = ((1./(disapearance_tmp./max(disapearance_tmp)))-1).*(sum(kernel.dispersal_wind_annulus(2:length(kernel.dispersal_wind_annulus))));
clear wind_edge_a
clear wind_edge_b
clear disapearance_tmp
Prop_among.wind_lost_within_source = 1-(sum(kernel.dispersal_wind_annulus) +  Prop_among.wind_lost + kernel.wind(1,1));
Prop_among.wind_export = sum(Prop_among.wind_annulus,1);

% Exclude wind dispersal if asked to.
if strcmpi('y',sensi_para.wind_dispersal_present) == 0
    Prop_among.wind_annulus(:,:) = 0;
    Prop_among.wind_edge(:) = 0;
    Prop_among.wind_lost = 0;
    Prop_among.wind_lost_within_source = 0;
end

%% calcute the human dispersal matrix throughout the meta_population from Kernel
for i = 1:pop_max
    for y = 1:1:pop_max
        x = Distance_matrix.centre(i,y);
        % we remove within patch dispersal from our calculations
        if x > 0
            Prop_among.human(i,y) = kernel.dispersal_human(find(Distance_matrix.oneD.centre == x));%#ok<*FNDSB>
        end
    end
end

% calculate two sided dissapearences
Prop_among.human_export = sum(Prop_among.human,1);
% exported is to two sites so for non edge patches it does double, this is
% correct since seeds that are dropped within the patch can have a
% tertiarty dispersal in another direction.
human_edge_1 = max(Prop_among.human_export)-Prop_among.human_export;
human_inf = fun_num((max(max(Distance_matrix.centre)))+step);
Prop_among.human_edge  = human_edge_1 + human_inf;
clear disapearance_tmp
clear x
clear i
clear y
% Exclude human dispersal if asked to.
if strcmpi('y',sensi_para.human_dispersal_present) == 0
    Prop_among.human(:,:) = 0;
    Prop_among.human_export(:) = 0;
    Prop_among.human_edge(:) = 0;
end

