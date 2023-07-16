function [Dc] = getGrad_contrainteT1etoile(W,spins,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

TR = getTR_changementvar(opt,W);
alpha = getalpha_changementvar(opt,W);


% c = opt.Kmin - cos(alpha)*exp(-TR/max(spins{1}.T1));

Dc = zeros(size(W));

if opt.optimisation_alpha
    Dc(opt.Np+1,1) = sin(alpha)*exp(-TR/max(spins{1}.T1));
end

if opt.optimisation_TR
    Dc(opt.Np+1,3) = 1/max(spins{1}.T1)*cos(alpha)*exp(-TR/max(spins{1}.T1));
end

%% Changement variable 
if opt.optimisation_TR
    if opt.changement_of_variable
        v = [W(1:end-1,3)' opt.Nlignes*W(end,3)];
        M1 = eye(opt.Np+1)*sum(v);
        M2 = diag(W(1:end,3)')*ones(opt.Np+1);
        M2(:,end)=opt.Nlignes*M2(:,end);
        D_ti_alphai = (M1-M2)/(sum(v)^2)*opt.tempsfixe_valeur;
    else
        D_ti_alphai = eye(opt.Np+1);
    end
else
    if opt.changement_of_variable
        v = [W(1:end,3)'];
        M1 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
        M1(1:opt.Np,1:opt.Np) = eye(opt.Np)*sum(v);
        M2 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
        M2(1:opt.Np,1:opt.Np) = diag(W(1:opt.Np,3)')*ones(opt.Np);
        D_ti_alphai = (M1-M2)/(sum(v)^2)*(opt.time_of_a_segment-opt.TR*opt.Nlignes);
    else
        D_ti_alphai = eye(size(W,1));
    end
end 

Dc(:,3) = D_ti_alphai'*Dc(:,3);

Dc = Dc(:);



end
