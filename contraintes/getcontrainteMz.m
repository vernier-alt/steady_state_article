function [contrainte] = getcontrainteMz(optimParam,spins,opt,f)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here


alpha = getalpha_changementvar(opt,optimParam);
TR = getTR_changementvar(opt,optimParam);
i = 1;
for p = 1:numel(opt.offsetVecHz):numel(spins) 
    
    if spins{p}.on_resonnance == false
        break;
    end
   
    EA = exp(-opt.TA/spins{p}.T1); %ms
    EB = exp(-opt.TB/spins{p}.T1); %ms

    E1 = exp(-TR/spins{p}.T1); %ms
    K = cos(alpha)*E1; %ms 

    H_U = [0 0 0 1]*spins{p}.U(:,:,opt.Np)*[1 0 0 0]';
    F_U =  [0 0 0 1]*spins{p}.U(:,:,opt.Np)*[0 0 0 1]';

    % STEADY - STATE 

        A = spins{p}.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
        lambda = K^opt.Nlignes*EA*EB;

        s = spins{p}.Mt0(end)*(1-E1)/(1-K);

        ss = ((A + lambda*H_U)/(1-lambda*F_U ));

        ss_EB = (ss -(1-EB))/EB;

        ff = opt.Nlignes - f +1; 
        
        if (spins{p}.saturation == true)

            contrainte(i) = (s + (ss_EB - s ) * K^(-ff))^2;
            i = i+1;

        end

end


end

