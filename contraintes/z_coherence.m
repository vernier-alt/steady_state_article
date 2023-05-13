function [contrainte] = z_coherence(optimParam,spins,opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

alpha = getalpha_changementvar(opt,optimParam);
TR = getTR_changementvar(opt,optimParam);
i = 1;

ez = [0 0 0 1]';

for p = 1:numel(opt.offsetVecHz):numel(spins) % on resonnance 
    
    if spins{p}.on_resonnance == false
        break;
    end
    
    A_0 = spins{p}.U(:,:,opt.Np);

    no1 = p+1 ;
    no2 = p+numel(opt.offsetVecHz)-1 ;
    
    for p_out = no1:no2 
        
            if spins{p_out}.on_resonnance == true
                break;
            end
            
            A_out = spins{p_out}.U(:,:,opt.Np);
%             contrainte(i) = (A_0*ez)'*(A_out*ez)/ ((A_0*ez)'*(A_0*ez))-0.8;
              contrainte(i) = 0.9*(A_0*ez)'*(A_0*ez)   - (A_0*ez)'*(A_out*ez);
            i = i+1;
    end
    
end


end

