function ti = getDelais_changementvar(opt,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if opt.changement_variable_delais
    ti = (opt.tempsfixe_valeur-opt.Nlignes*opt.TR)*W(1:opt.Np,3)/(sum(W(1:opt.Np,3)));
else
    ti = W(:,3);
end

end
