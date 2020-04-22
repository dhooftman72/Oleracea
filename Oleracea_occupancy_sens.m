if strcmpi('y',Do_sensitivity) ~= 0
    mini = Output.Sensi_para_baseline;
else

    for sensis = 1:1:length(Sensis)
        i = Sensis(sensis);
        max_length = length(Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean);
        test(sensis) =  (Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean(max_length))-...
            (Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean(1));
    end
        limit = 0.5.*Output.s1.Occupancy.adult.mean(1);
    test(test<limit) = max(test);
    mini = Sensis(find(test == min(test)));
    mini= mini(1);
end


baseline_occu = min(find((Output.(genvarname(['s',int2str(mini)])).Occupancy.adult.mean) >...
(Output.(genvarname(['s',int2str(mini)])).Occupancy.adult.mean(1).*1.5)));
 if isempty(baseline_occu) == 1
            baseline_occu = NaN;
 end
for sensis = 1:1:length(Sensis)
    i = Sensis(sensis);
    occu_doubles= min(find(Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean >...
        (Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean(1)*1.5)));
    if isempty( occu_doubles) == 1
        occu_doubles= min(find(Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean <...
            (Output.(genvarname(['s',int2str(i)])).Occupancy.adult.mean(1)*(1/1.5))));
        if isempty(occu_doubles) == 1
            occu_doubles = NaN;
        end
    end
    base = baseline_occu;
    new = occu_doubles;
    if isnan(base) ~= 1 && isnan(new) ~= 1 % This is SENSITIVITY !!!
        %Lambda_meta   
        sensitivity = base/new; %(base-new)/base;
        if occu_doubles == 100;
            sensitivity = 1;
        end
%         if new < base
%             sensitivity = new/base;
%         end

        %    sensitivity = sensitivity./abs(1-change_factor);
        Output.Occupancy_meta(sensis,10) = sensitivity;
        Output.Occupancy_meta(sensis,11) =  occu_doubles;
        Output.Occupancy_meta(sensis,12) = Output.(genvarname(['s',int2str(Sensis(sensis))])).Occupancy.adult.std(new) ;
    end
    
end