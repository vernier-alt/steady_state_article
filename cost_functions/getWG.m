function [WGSD] = getWG(optimParam,spins,opt)
%UNTITLED9 Summary of this function goes here
%  C = - |Mmax| + | Mmin| 

WGSD = 0;

for p = 1:numel(spins) 
    
    EA = exp(-opt.TA/spins{p}.T1); %ms
    EB = exp(-opt.TB/spins{p}.T1); %ms
    
    TR = opt.tempsfixe_valeur*optimParam(end,3)/(sum(optimParam(1:end-1,3))+opt.Nlignes*optimParam(end,3));
    E1 = exp(-TR/spins{p}.T1); %ms
    
    K = cos(optimParam(end,1))*E1; %ms 

    H_U = [0 0 0 1]*spins{p}.U(:,:,opt.Np)*[1 0 0 0]';
    F_U =  [0 0 0 1]*spins{p}.U(:,:,opt.Np)*[0 0 0 1]';

    % STEADY - STATE 

    A = spins.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
    lambda = K^opt.Nlignes*EA*EB;

    s = spins.Mt0(end)*(1-E1)/(1-K);

    ss = ((A + lambda*H_U)/(1-lambda*F_U ));

    ss_EB = (ss -(1-EB))/EB;

    for j =1:numel(opt.vec)
        
        f = opt.Nlignes - opt.vec(j); 
        
        if (spins{p}.max == -1)
            
            WGSD = WGSD +abs((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(optimParam(end,1)));

        else     
            WGSD = WGSD -abs((s + (ss_EB - s ) * K^(-f))*exp(-opt.TE/spins{p}.T2)*sin(optimParam(end,1)));

        end

    end
end

end