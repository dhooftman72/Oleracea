if sensi == 2 % all SB
     SBincorp = SBincorp .*sensi_para.change_factor;
     SBrate = SBrate.*sensi_para.change_factor;
        % function needs to be within limits 0-100%
        if SBincorp < 0
            SBincorp = 0;
        end
        if SBincorp > 1
            SBincorp = 1;
        end
        if SBrate < 0
           SBrate = 0;
        end
        if SBrate> 1
            SBrate = 1;
        end
elseif sensi == 3 % All GS
    for f= 1:1:length(Gs_dis)
        % function needs to be within limits 0-100%
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
    end
elseif sensi == 4
    for f= 1:1:length(Ss_dis)
        % function needs to be within limits 0-100%
        Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end
    end
elseif sensi == 5    
     for i = 1:1:max_diameter
        Adult_surv(2,i) = Adult_surv(2,i).*sensi_para.change_factor;      
        % function needs to be within limits 0-100%
        if  Adult_surv(2,i) < 0
            Adult_surv(2,i) = 0;
        end
        if  Adult_surv(2,i) >1
            Adult_surv(2,i) = 1;
        end
     end
elseif sensi == 6        
     for i = 1:1:max_diameter
        Adult_increase(2,i) = Adult_increase(2,i).*sensi_para.change_factor;  
     end
     % function NOT needs to be within limits 0-100%
elseif sensi == 7    
     for i = 1:1:max_diameter
        Flowering_like(2,i) = Flowering_like(2,i).*sensi_para.change_factor;      
        % function needs to be within limits 0-100%
        if  Flowering_like(2,i) < 0
            Flowering_like(2,i) = 0;
        end
        if  Flowering_like(2,i) >1
            Flowering_like(2,i) = 1;
        end
     end   
elseif sensi == 8        
    for i = 1:1:max_diameter
        Fecundity(2,i) = Fecundity(2,i).*sensi_para.change_factor;  
     end
     % function NOT needs to be within limits 0-100%
elseif sensi == 9
    sensi_para.Pu = sensi_para.change_factor;
elseif sensi == 10
   sensi_para.Pu = sensi_para.change_factor;
elseif sensi == 11
% is done in define file
elseif sensi == 12
% is done in define file
elseif sensi == 13 % all SB
     SBincorp = SBincorp .*sensi_para.change_factor;
        % function needs to be within limits 0-100%
        if SBincorp < 0
            SBincorp = 0;
        end
        if SBincorp > 1
            SBincorp = 1;
        end
elseif sensi == 14 % all SB
     SBrate = SBrate.*sensi_para.change_factor;
        % function needs to be within limits 0-100%
        if SBrate < 0
           SBrate = 0;
        end
        if SBrate> 1
            SBrate = 1;
        end
elseif sensi == 15 % All GS
        % function needs to be within limits 0-100%
        f = 1;
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
elseif sensi == 16 % All GS
        % function needs to be within limits 0-100%
        f = 2;
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
elseif sensi == 17 % All GS
        % function needs to be within limits 0-100%
        f = 3;
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
elseif sensi == 18 % All GS
        % function needs to be within limits 0-100%
        f = 4;
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
elseif sensi == 19 % All GS
        % function needs to be within limits 0-100%
        f = 5;
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
elseif sensi== 20
    f = 1;    
    Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end
elseif sensi== 21
    f = 2;    
    Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end
elseif sensi== 22
    f = 3;    
    Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end
elseif sensi== 23
    f = 4;    
    Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end     
elseif sensi== 24
    f = 5;    
    Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end   
elseif sensi == 25 % GS 5 en 3 combined
        % function needs to be within limits 0-100%
        Gs_dis(3) = Gs_dis(3).*sensi_para.change_factor;
        if Gs_dis(3) < 0
            Gs_dis(3) = 0;
        end
        if Gs_dis(3) > 1
            Gs_dis(3) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(3) = Gs_undis(3).*sensi_para.change_factor;
        if Gs_undis(3) < 0
            Gs_undis(3) = 0;
        end
        if Gs_undis(3) > 1
            Gs_dis(3) = 1;
        end
        
        Gs_dis(5) = Gs_dis(5).*sensi_para.change_factor;
        if Gs_dis(5) < 0
            Gs_dis(5) = 0;
        end
        if Gs_dis(5) > 1
            Gs_dis(5) = 1;
        end
        % function needs to be within limits 0-100%
        Gs_undis(5) = Gs_undis(5).*sensi_para.change_factor;
        if Gs_undis(5) < 0
            Gs_undis(5) = 0;
        end
        if Gs_undis(5) > 1
            Gs_dis(5) = 1;
        end
elseif sensi == 26    
     for i = 1:1:9
        Adult_surv(2,i) = Adult_surv(2,i).*sensi_para.change_factor;      
        % function needs to be within limits 0-100%
        if  Adult_surv(2,i) < 0
            Adult_surv(2,i) = 0;
        end
        if  Adult_surv(2,i) >1
            Adult_surv(2,i) = 1;
        end
     end        
elseif sensi == 27    
     for i = 10:1:max_diameter
        Adult_surv(2,i) = Adult_surv(2,i).*sensi_para.change_factor;      
        % function needs to be within limits 0-100%
        if  Adult_surv(2,i) < 0
            Adult_surv(2,i) = 0;
        end
        if  Adult_surv(2,i) >1
            Adult_surv(2,i) = 1;
        end
     end 
elseif sensi == 28 % All GS
    for f= 1:1:length(Gs_dis)
        % function needs to be within limits 0-100%
        Gs_dis(f) = Gs_dis(f).*sensi_para.change_factor;
        if Gs_dis(f) < 0
            Gs_dis(f) = 0;
        end
        if Gs_dis(f) > 1
            Gs_dis(f) = 1;
        end
    end
elseif sensi == 29 % All GS
    for f= 1:1:length(Gs_undis)
        % function needs to be within limits 0-100%
        Gs_undis(f) = Gs_undis(f).*sensi_para.change_factor;
        if Gs_undis(f) < 0
            Gs_undis(f) = 0;
        end
        if Gs_undis(f) > 1
            Gs_dis(f) = 1;
        end
    end
elseif sensi == 30
    for f= 1:1:length(Ss_dis)
        % function needs to be within limits 0-100%
        Ss_dis(f) = Ss_dis(f).*sensi_para.change_factor;
        if Ss_dis(f) < 0
            Ss_dis(f) = 0;
        end
        if Ss_dis(f) > 1
            Ss_dis(f) = 1;
        end
    end
elseif sensi == 31
    for f= 1:1:length(Ss_undis)        
        Ss_undis(f) = Ss_undis(f).*sensi_para.change_factor;
        if Ss_undis(f) < 0
            Ss_undis(f) = 0;
        end
        if Ss_undis(f) > 1
            Ss_undis(f) = 1;
        end
    end
elseif sensi == 32
        % function needs to be within limits 0-100%
        Ss_dis(2) = Ss_dis(2).*sensi_para.change_factor;
        if Ss_dis(2) < 0
            Ss_dis(2) = 0;
        end
        if Ss_dis(2) > 1
            Ss_dis(2) = 1;
        end
        % function needs to be within limits 0-100%
        Ss_undis(2) = Ss_undis(2).*sensi_para.change_factor;
        if Ss_undis(2) < 0
            Ss_undis(2) = 0;
        end
        if Ss_undis(2) > 1
            Ss_dis(2) = 1;
        end
        
        Ss_dis(3) = Ss_dis(3).*sensi_para.change_factor;
        if Ss_dis(3) < 0
            Ss_dis(3) = 0;
        end
        if Ss_dis(3) > 1
            Ss_dis(3) = 1;
        end
        % function needs to be within limits 0-100%
        Ss_undis(3) = Ss_undis(3).*sensi_para.change_factor;
        if Ss_undis(3) < 0
            Ss_undis(3) = 0;
        end
        if Ss_undis(3) > 1
            Ss_dis(3) = 1;
        end
        
        Ss_dis(4) = Ss_dis(4).*sensi_para.change_factor;
        if Ss_dis(4) < 0
            Ss_dis(4) = 0;
        end
        if Ss_dis(4) > 1
            Ss_dis(4) = 1;
        end
        % function needs to be within limits 0-100%
        Ss_undis(4) = Ss_undis(4).*sensi_para.change_factor;
        if Ss_undis(4) < 0
            Ss_undis(4) = 0;
        end
        if Ss_undis(4) > 1
            Ss_dis(4) = 1;
        end
        
        Ss_dis(5) = Ss_dis(5).*sensi_para.change_factor;
        if Ss_dis(5) < 0
            Ss_dis(5) = 0;
        end
        if Ss_dis(5) > 1
            Ss_dis(5) = 1;
        end
        % function needs to be within limits 0-100%
        Ss_undis(5) = Ss_undis(5).*sensi_para.change_factor;
        if Ss_undis(5) < 0
            Ss_undis(5) = 0;
        end
        if Ss_undis(5) > 1
            Ss_dis(5) = 1;
        end
        
        
end
    