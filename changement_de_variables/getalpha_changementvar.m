function alpha = getalpha_changementvar(opt,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if opt.optimisation_alpha
    alpha = W(end,1);
else

    alpha = opt.alpha;
end 

end

