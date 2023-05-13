function M = getRelaxMat(spins,tp,expF,mode)
 
if strcmp(mode,'exp')
    
B0 = spins.B0_inh ;

M  = expF([0 0 0 0 ; ...
    0 -1/spins.T2 B0 0 ; ...
    0 -B0 -1/spins.T2 0 ; ...
    spins.Mt0(end)/spins.T1 0 0 -1/spins.T1]...
    *tp);

% changement sens
% M  = expF([0                    0               0               0 ; ...
%            0                    -1/spins.T2     -B0              0 ; ...
%            0                    B0             -1/spins.T2     0 ; ...
%     spins.Mt0(end)/spins.T1     0               0       -1/spins.T1]...
%     *tp);




elseif strcmp(mode,'exact')
    

    c = cos(spins.B0_inh*tp);
    s = sin(spins.B0_inh*tp);
  
    M = [1                                  0                           0                       0;...
         0                                  exp(-tp/spins.T2)*c       -exp(-tp/spins.T2)*s      0;...
         0                                  exp(-tp/spins.T2)*s        exp(-tp/spins.T2)*c      0;...
         spins.Mt0(end)*(1-exp(-tp/spins.T1))   0                           0                       exp(-tp/spins.T1)];
      


end

end