function [signal]= get_signal(optimParam,spins,opt,ll)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here 

% EXP
f = opt.Nlignes - ll+1;
% disp('get_signal');
    EA = exp(-opt.TA/spins.T1); %ms
    EB = exp(-opt.TB/spins.T1); %ms
    TR = opt.TR;
    alpha = opt.alpha;
    E1 = exp(-TR/spins.T1); %ms
    K = cos(alpha)*E1; %ms 

    H_U = [0 0 0 1]*spins.U(:,:,opt.Np)*[1 0 0 0]';
    F_U =  [0 0 0 1]*spins.U(:,:,opt.Np)*[0 0 0 1]';

    % STEADY - STATE 

        A = spins.Mt0(end)* ((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
        lambda = K^opt.Nlignes*EA*EB;

        s = spins.Mt0(end)*(1-E1)/(1-K);

        ss = ((A + lambda*H_U)/(1-lambda*F_U ));

        ss_EB = (ss -(1-EB))/EB;
        
        signal = ((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins.T2)*sin(alpha));
            
end


















