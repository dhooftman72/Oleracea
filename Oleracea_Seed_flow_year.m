Seed_present = Seeds_total;
val = rand(pop_max,4);
%% first wind dipersal
for population = 1:pop_max
        % lost
        Seed_flow_tmp = Seed_present(population).*Prop_among.wind_lost;
        Seed_flow_tmp(Seed_flow_tmp<val(population,1) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Lost_seeds(population)= sum(Seed_flow_tmp);
        clear Seed_flow_tmp
        
        Seed_flow_tmp = Seed_present(population).*Prop_among.wind_lost_within_source;
        Seed_flow_tmp(Seed_flow_tmp<val(population,2) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Lost_within_seeds(population)= sum(Seed_flow_tmp);
        clear Seed_flow_tmp
       
        Seed_flow_tmp = Seed_present(population).*Prop_among.wind_edge(population);
        Seed_flow_tmp(Seed_flow_tmp<val(population,3) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Edge_seeds(population)= sum(Seed_flow_tmp);
        clear Seed_flow_tmp
        
        Seed_flow_tmp = Seed_present(population).*Prop_among.wind_export(population);
        Seed_flow_tmp(Seed_flow_tmp<val(population,4) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Export_seeds(population)= sum(Seed_flow_tmp);
        
        % Import 
        Seed_flow_tmp = Seed_present(population).*(Prop_among.wind_annulus(population,:));
        Seed_flow_tmp(Seed_flow_tmp<val(population) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Recieved_seeds_tmp(population,:)= Seed_flow_tmp;
        clear Seed_flow_tmp
end
Recieved_seeds = sum(Recieved_seeds_tmp,1)
Exported_seeds = (Export_seeds + Lost_seeds + Edge_seeds + Lost_within_seeds)
Seed_present_w = Seed_present - Exported_seeds + Recieved_seeds;
clear val

%% human dispersal
% calculate the number of seeds that will be available on the patch and
% picked up.  
prop_area = sensi_para.area_size/((sensi_para.area_size)^2); % note the area is considered square and the path 1 m wide  
total_available = Seed_present_w.* prop_area.* (sensi_para.Pu.* sensi_para.pick_up_sensitivity);

val = rand(size(Prop_among.human,1),1);
for population = 1:pop_max
        Seed_flow_tmp = total_available(population)'.*(Prop_among.human(population,:));
        Seed_flow_tmp(Seed_flow_tmp<val(population) & Seed_flow_tmp<1) = 0;
        Seed_flow_tmp(Seed_flow_tmp>0 & Seed_flow_tmp<1) = 1;
        Seed_flow_h(population,:)= Seed_flow_tmp;
        
        Seed_flow_export(population) =  Prop_among.human_export(population).*total_available(population);
        Seed_flow_edge(population) = Prop_among.human_edge(population).*total_available(population);
        clear test
        clear Seed_flow_tmp
end
clear val
%Recalculate the seeds after dispersal
Recieved_seeds_h = sum(Seed_flow_h,1);
Exported_seeds_h = Seed_flow_export + Seed_flow_edge;
Seed_present_h = Seed_present_w - Exported_seeds_h + Recieved_seeds_h;
Seeds_total = Seed_present_h;
clear Seed_present_h
clear Seed_present_w 
clear Exported_seeds_h
clear Recieved_seeds_h
clear Exported_seeds
clear Recieved_seeds
clear Seed_flow
clear Seed_flow_h
clear Seed_present
clear Seed_total_old