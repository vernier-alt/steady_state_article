function ti = getDelais_changementvar(opt,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if opt.changement_of_variable
    ti = (opt.time_of_a_segment-opt.Nlignes*opt.TR)*W(1:opt.Np,3)/(sum(W(1:opt.Np,3)));
else
    ti = W(:,3);
end

end
