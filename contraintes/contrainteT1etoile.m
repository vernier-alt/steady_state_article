function [c] = contrainteT1etoile(W,spins,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

TR = getTR_changementvar(opt,W);
alpha = getalpha_changementvar(opt,W);

c = opt.Kmin - cos(alpha)*exp(-TR/max(spins{1}.T1));

end

