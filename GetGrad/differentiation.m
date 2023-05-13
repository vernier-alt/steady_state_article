function [D_final] = differentiation(W,spins,opt,option)
% AVEC ti

TE = opt.TE;
TA = opt.TA;
TB = opt.TB;
T1 = spins.T1; 
T2 = spins.T2; 
n=opt.Nlignes;
M0 = spins.Mt0(end);
% coeff = opt.coeff;

TR = getTR_changementvar(opt,W);
alpha = getalpha_changementvar(opt,W);

E1 = exp(-TR/T1);
EB = exp(-TB/T1);
EA = exp(-TA/T1);
K = cos(alpha)*E1; %ms

H_U = [0 0 0 1]*spins.U(:,:,opt.Np)*[1 0 0 0]';
F_U = [0 0 0 1]*spins.U(:,:,opt.Np)*[0 0 0 1]';

A = M0*((1-EB)+( (1-E1)*(1-K^n)/(1-K)+K^n*(1-EA))*EB);
lambda = K^n*EA*EB;
s = M0*(1-E1)/(1-K);
sz = ((A + lambda*H_U)/(1-lambda*F_U ));
ss_EB = (sz -(1-EB))/EB;

D_final = zeros(3,1);

for jj=1:numel(opt.vec)

    f = opt.Nlignes - opt.vec(jj)+1;
    signal = (s + (ss_EB - s ) * K^(-f))*exp(-TE/T2)*sin(alpha);
    Mz = (s + (ss_EB - s ) * K^(-f));
    
    D_Mz_alpha = (M0*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2 - (exp(TB/T1)*((M0*exp(-TB/T1)*(n*exp(-TR/T1)*sin(alpha)*(exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1) - (exp(-TR/T1)*sin(alpha)*((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2 + (n*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - H_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - (F_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*sin(alpha)*(M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1)^2) + (M0*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2)/(exp(-TR/T1)*cos(alpha))^f + (f*exp(-TR/T1)*sin(alpha)*(exp(TB/T1)*(exp(-TB/T1) + (M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - 1) - (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)))/(exp(-TR/T1)*cos(alpha))^(f + 1);
    D_Mz_TR = (f*exp(-TR/T1)*cos(alpha)*(exp(TB/T1)*(exp(-TB/T1) + (M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - 1) - (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)))/(T1*(exp(-TR/T1)*cos(alpha))^(f + 1)) - (M0*exp(-TR/T1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) - (exp(TB/T1)*((M0*exp(-TB/T1)*((exp(-TR/T1)*((exp(-TR/T1)*cos(alpha))^n - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) + (n*exp(-TR/T1)*cos(alpha)*(exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/T1 - (exp(-TR/T1)*cos(alpha)*((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2) + (n*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1))) - (H_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1)*cos(alpha))^(n - 1))/T1)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - (F_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*cos(alpha)*(M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(T1*(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1)^2)) - (M0*exp(-TR/T1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) + (M0*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2))/(exp(-TR/T1)*cos(alpha))^f + (M0*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2);
 

    
    
    D_alpha = exp(-TE/T2)*cos(alpha)*((exp(TB/T1)*(exp(-TB/T1) + (M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - 1) - (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1))/(exp(-TR/T1)*cos(alpha))^f + (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) + exp(-TE/T2)*sin(alpha)*((M0*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2 - (exp(TB/T1)*((M0*exp(-TB/T1)*(n*exp(-TR/T1)*sin(alpha)*(exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1) - (exp(-TR/T1)*sin(alpha)*((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2 + (n*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - H_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - (F_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*sin(alpha)*(M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1)^2) + (M0*exp(-TR/T1)*sin(alpha)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)^2)/(exp(-TR/T1)*cos(alpha))^f + (f*exp(-TR/T1)*sin(alpha)*(exp(TB/T1)*(exp(-TB/T1) + (M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - 1) - (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)))/(exp(-TR/T1)*cos(alpha))^(f + 1));
    D_TR = -exp(-TE/T2)*sin(alpha)*((exp(TB/T1)*((M0*exp(-TB/T1)*((exp(-TR/T1)*((exp(-TR/T1)*cos(alpha))^n - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) + (n*exp(-TR/T1)*cos(alpha)*(exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/T1 - (exp(-TR/T1)*cos(alpha)*((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2) + (n*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1))) - (H_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1)*cos(alpha))^(n - 1))/T1)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - (F_U*n*exp(-TA/T1)*exp(-TB/T1)*exp(-TR/T1)*cos(alpha)*(M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)*(exp(-TR/T1)*cos(alpha))^(n - 1))/(T1*(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1)^2)) - (M0*exp(-TR/T1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) + (M0*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2))/(exp(-TR/T1)*cos(alpha))^f + (M0*exp(-TR/T1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)) - (f*exp(-TR/T1)*cos(alpha)*(exp(TB/T1)*(exp(-TB/T1) + (M0*(exp(-TB/T1) + exp(-TB/T1)*((exp(-TA/T1) - 1)*(exp(-TR/T1)*cos(alpha))^n + (((exp(-TR/T1)*cos(alpha))^n - 1)*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)) - 1) - H_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n)/(F_U*exp(-TA/T1)*exp(-TB/T1)*(exp(-TR/T1)*cos(alpha))^n - 1) - 1) - (M0*(exp(-TR/T1) - 1))/(exp(-TR/T1)*cos(alpha) - 1)))/(T1*(exp(-TR/T1)*cos(alpha))^(f + 1)) - (M0*exp(-TR/T1)*cos(alpha)*(exp(-TR/T1) - 1))/(T1*(exp(-TR/T1)*cos(alpha) - 1)^2));
    
    if strcmp(option,'contrainte')
     
        if opt.optimisation_alpha
            D_final(1,1) = D_final(1,1) + 2*Mz*D_Mz_alpha;
        end
        if opt.optimisation_TR
            D_final(3,1) = D_final(3,1) + 2*Mz*D_Mz_TR ;  
        end
     
    elseif strcmp('@getWG_epsilon',opt.costFunction)

      D_norme = signal*power(opt.epsilon_abs^2 + signal.^2,-0.5);


            if (spins.max == -1)&& (spins.saturation == false)
                if opt.optimisation_alpha
                    D_final(1,1) = D_final(1,1) + D_alpha*D_norme;
                end
                if opt.optimisation_TR
                    D_final(3,1) = D_final(3,1) + D_TR*D_norme;  
                end
            else 
                if opt.optimisation_alpha
                    D_final(1,1) = D_final(1,1) +  D_alpha*D_norme;
                end
                if opt.optimisation_TR
                    D_final(3,1) = D_final(3,1) + D_TR*D_norme;
                    
                end
            end
            
            
            
     end 
      
     %       if strcmp('@getCO',opt.costFunction) 
% 
%             if (spins.max == -1)&& (spins.saturation == false)
%                 if opt.optimisation_alpha
%                     D_final(1,1) = D_final(1,1) + 2*opt.sat*signal*D_alpha;
%                 end
%                 if opt.optimisation_TR
%                     D_final(3,1) = D_final(3,1) + 2*opt.sat*signal*D_TR;  
%                 end
%             else 
%                 if opt.optimisation_alpha
%                     D_final(1,1) = D_final(1,1) +  2*signal*D_alpha;
%                 end
%                  if opt.optimisation_TR
%                     D_final(3,1) = D_final(3,1) + 2*signal*D_TR;
%                     
%                  end
%             end
% 
%       elseif strcmp('@getWGABS',opt.costFunction)
% 
%             D_final(1,1) = D_final(1,1) + sign(signal)*D_alpha;
%             D_final(3,1) = D_final(3,1) + sign(signal)*D_TR; 

end


