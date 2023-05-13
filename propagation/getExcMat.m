function M = getExcMat(opt,spins,optimparam)

if strcmp(opt.mode,'exp')

              
    EXPMAT  = str2func( opt.expmFunction ) ;

    M = EXPMAT( [0                      0                           0                           0 ; ...
                 0                      -1/spins.T2                 spins.B0_inh                -optimparam(2)*spins.B1_inh ; ...
                 0                      -spins.B0_inh               -1/spins.T2                 optimparam(1)*spins.B1_inh ;  ...
                 spins.Mt0(end)/spins.T1 optimparam(2)*spins.B1_inh -optimparam(1)*spins.B1_inh -1/spins.T1 ]...
                  *opt.tp ) ;        
              
%     M = EXPMAT( [0                          0                               0                                   0 ; ...
%                  0                          -1/T2                     -Omegaz                       optimparam(2)*spins.B1_inh ; ...
%                  0                          Omegaz                    -1/T2                          -optimparam(1)*spins.B1_inh ;  ...
%                  spins.Mt0(end)/T1    -optimparam(2)*spins.B1_inh    optimparam(1)*spins.B1_inh         -1/T1 ]...
%                   *opt.tp ) ;






elseif strcmp(opt.mode,'exact')

%     if Omegaz ==0
        
%         b = optimparam(1); % positif
%         c = cos(-b*opt.tp);
%         s = sin(-b*opt.tp);
% 
%         optimparamx_n = cos(optimparam(2));
%         optimparamy_n = sin(optimparam(2));
% 
%         M =      [1     0                                                         0                                                   0
%                   0     (optimparamx_n)^2*(1-c)+c                         (optimparamx_n*optimparamy_n)*(1-c)                        optimparamy_n*s;...
%                   0     optimparamx_n*optimparamy_n*(1-c)                         (optimparamy_n)^2*(1-c)+c                          -optimparamx_n*s;
%                   0     -optimparamy_n*s                                   optimparamx_n*s                                       c   ];

        b = optimparam(1); % positif
        c = cos(b*opt.tp);
        s = sin(b*opt.tp);

        optimparamx_n = cos(optimparam(2));
        optimparamy_n = sin(optimparam(2));

        M =      [1     0                                                         0                                                   0
                  0     (optimparamx_n)^2*(1-c)+c                         (optimparamx_n*optimparamy_n)*(1-c)                        optimparamy_n*s;...
                  0     optimparamx_n*optimparamy_n*(1-c)                         (optimparamy_n)^2*(1-c)+c                          -optimparamx_n*s;
                  0     -optimparamy_n*s                                   optimparamx_n*s                                       c   ];


    
%     end
end

end


% 
% M = EXPMAT( [0                           0                       0                      0           ; ...

%              0                       -1/T2            Omegaz        -optimparam(2)*spins.optimparam_inh ; ...

%              0                      -Omegaz           -1/T2          optimparam(1)*spins.optimparam_inh ;  ...  

%     spins.Mt0(end)/T1       optimparam(2)*spins.optimparam_inh       -optimparam(1)*spins.optimparam_inh     -1/T1      ]  *opt.tp ) ;
%      


