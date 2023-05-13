function [ss] = get_ss(optimParam,spins,opt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here


% EXP

% Plusieurs espèces 
ss = 0;
for s = 1:numel(spins) 

    EA = exp(-opt.TA/spins{s}.T1); %ms
    EB = exp(-opt.TB/spins{s}.T1); %ms
    TR = opt.tempsfixe_valeur*optimParam(end,3)/(sum(optimParam(1:end-1,3))+opt.Nlignes*optimParam(end,3));
    E1 = exp(-TR/spins{s}.T1); %ms
    K = cos(optimParam(end,1))*E1; %ms 

    H_U = [0 0 0 1]*spins{s}.U(:,:,opt.Np)*[1 0 0 0]';
    F_U =  [0 0 0 1]*spins{s}.U(:,:,opt.Np)*[0 0 0 1]';

    % STEADY - STATE 

    A = spins.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
    lambda = K^opt.Nlignes*EA*EB;


    ss = ss + ((A + lambda*H_U)/(1-lambda*F_U ));

end

end

