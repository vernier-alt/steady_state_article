function ti = getDelais_changementvar(opt,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if opt.optimisation_TR
    if opt.changement_variable_delais
        ti = opt.tempsfixe_valeur*W(:,3)/(sum(W(1:end-1,3))+opt.Nlignes*W(end,3));
    else
        ti = W(:,3);
    end
else
    if opt.changement_variable_delais
        ti = (opt.tempsfixe_valeur-opt.Nlignes*opt.TR)*W(1:opt.Np,3)/(sum(W(1:opt.Np,3)));
    else
        ti = W(:,3);
    end

end 

end
