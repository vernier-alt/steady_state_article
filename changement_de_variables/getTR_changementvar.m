function TR = getTR_changementvar(opt,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if opt.optimisation_TR
    if opt.changement_variable_delais
        TR = opt.tempsfixe_valeur*W(end,3)/(sum(W(1:end-1,3))+opt.Nlignes*W(end,3));
    else
        TR = W(end,3);
    end
else % ce n'est pas une variable d'optimisation
    TR = opt.TR;
end


end

