function [WGSD] = getCO(optimParam,spins,opt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% EXP

numerateur = 0;

alpha = getalpha_changementvar(opt,optimParam);
TR = getTR_changementvar(opt,optimParam);

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
                
                numerateur = numerateur + opt.sat*((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))^2;
            
            elseif    (spins{p}.saturation == false)
%                 numerateur = numerateur -sign((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))*((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))^2;
                numerateur = numerateur -((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(alpha))^2;
                
            end
        end
end

WGSD = numerateur;

end

