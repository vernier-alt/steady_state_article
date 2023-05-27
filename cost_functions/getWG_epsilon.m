function [WGSD] = getWG_epsilon(optimParam,spins,opt)
% Cost function
%  C = - SQRT( Mmax^2 + e^2) + opt.sat*SQRT( Mmin^2 + e^2) or - SQRT( Mmax^2 + e^2) if constraint Mmin = 0 at the echo at line(s) opt.vec 

WGSD = 0;

alpha = opt.alpha;
TR = opt.TR;

for p = 1:numel(spins) 

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
    
    for j =1:numel(opt.vec)

        f = opt.Nlignes - opt.vec(j)+1;        
        if (spins{p}.max == -1) && (spins{p}.saturation == false)

            WGSD = WGSD + opt.sat*sqrt(((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))^2 + opt.epsilon_abs^2);

        elseif    (spins{p}.saturation == false)

            WGSD = WGSD - sqrt(((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))^2 + opt.epsilon_abs^2);

        end

    end

end

end

