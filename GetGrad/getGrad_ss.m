function g = getGrad_ss(W,spins,opt)

% Compute gradient dC/dw with complex differentiation (see M. Lapert paper)

EXPMAT  = str2func(opt.expmFunction) ;   

ez = [0 0 0 1]';
e0 = [1 0 0 0]';

% sum = 0;

    if opt.changement_variable_delais
        
        v = [W(1:end,3)'];
        M1 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
        M1(1:opt.Np,1:opt.Np) = eye(opt.Np)*sum(v);
        M2 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
        M2(1:opt.Np,1:opt.Np) = diag(W(1:opt.Np,3)')*ones(opt.Np);
        D_ti_alphai = (M1-M2)/(sum(v)^2)*(opt.tempsfixe_valeur-opt.TR*opt.Nlignes);

    else
        D_ti_alphai = eye(size(W,1));
    end
alpha = opt.alpha;
TR = opt.TR;

% cette fction de cout necessite la somme des spins d'un même tissue pour
% differentes valeurs de w0
somme = zeros(opt.N_iso,1);
if strcmp('@getWG_epsilon_B0',opt.costFunction)
    for tissue = 1:opt.N_iso % differents tissus
        for W0 = 1:numel(opt.offsetVecHz)
            num = (tissue - 1)*numel(opt.offsetVecHz) + W0;
            if (spins{num}.saturation == false) % pas besoin de calculer pour les spins a saturer
                EA = exp(-opt.TA/spins{num}.T1); %ms
                EB = exp(-opt.TB/spins{num}.T1); %ms
                E1 = exp(-TR/spins{num}.T1); %ms
                K = cos(alpha)*E1; %ms 
                A = spins{num}.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
                lambda = K^opt.Nlignes*EA*EB;
                s = spins{num}.Mt0(end)*(1-E1)/(1-K);
                H_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*e0;
                F_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*[0 0 0 1]';
                somme(tissue) = somme(tissue) + (s + ( ((A + lambda*H_U)/(1-lambda*F_U)*(EB^-1)-(EB^-1 - 1))-s)*K^(-(opt.Nlignes-opt.vec+1)))*(exp(-opt.TE/spins{num}.T2)*sin(alpha));
            end
        end
    end
end

% for num =1:numel(opt.offsetVecHz):numel(spins)
% for num =1:numel(spins)
for tissue = 1:opt.N_iso % differents tissus
    for W0 = 1:numel(opt.offsetVecHz)
        
        num = (tissue - 1)*numel(opt.offsetVecHz) + W0;
        
        ez_it = [0 0 0 1]';
        EA = exp(-opt.TA/spins{num}.T1); %ms
        EB = exp(-opt.TB/spins{num}.T1); %ms
        E1 = exp(-TR/spins{num}.T1); %ms
        K = cos(alpha)*E1; %ms 
        A = spins{num}.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
        lambda = K^opt.Nlignes*EA*EB;
        H_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*e0;
        F_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*[0 0 0 1]';
        s = spins{num}.Mt0(end)*(1-E1)/(1-K);
        signal = (s + ( ((A + lambda*H_U)/(1-lambda*F_U)*(EB^-1)-(EB^-1 - 1))-s)*K^(-(opt.Nlignes-opt.vec+1)))*(exp(-opt.TE/spins{num}.T2)*sin(alpha));
        fact = (EB^-1)*K^(-(opt.Nlignes-opt.vec+1))*(exp(-opt.TE/spins{num}.T2)*sin(alpha));

        for p = opt.Np:-1:1


            if strcmp(opt.mode,'exp') % B1xBy

                optParam_b1   = W(p,:); 
                del = getDelais_changementvar(opt,W);
                optParam_b1(3) = del(p);
                optParam_b1(1:2) = optParam_b1(1:2)*opt.mu;

                optParam_b1x = optParam_b1 + [1i*opt.epsilon 0 0];
                optParam_b1y = optParam_b1 + [0 1i*opt.epsilon 0];
                optParam_De     = optParam_b1 + [0 0 1i*opt.epsilon];     


                 % Relax_mat =getRelaxMat( spins{s},optParam_b1x(1,3),EXPMAT ) ;
                mat_iB1x = getRelaxMat( spins{num},optParam_b1x(1,3),EXPMAT ,opt.mode) * getExcMat(opt,spins{num},optParam_b1x);
                mat_iB1y = getRelaxMat( spins{num},optParam_b1y(1,3),EXPMAT ,opt.mode) * getExcMat(opt,spins{num},optParam_b1y);
                mat_iDe  = getRelaxMat( spins{num},optParam_De(1,3),EXPMAT ,opt.mode) * getExcMat(opt,spins{num},optParam_De );

                ez_b1x    = ez_it'*mat_iB1x ;
                ez_b1y    = ez_it'*mat_iB1y ;
                ez_De     = ez_it'*mat_iDe ;

                ez_it   = (real(ez_b1x + ez_b1y + ez_De)/3)' ; % 


                if p>1
                    spins{num}.gh(p,1) = 1/opt.epsilon*imag(ez_b1x*spins{num}.U(:,:,p-1)*e0) ; 
                    spins{num}.gh(p,2) = 1/opt.epsilon*imag(ez_b1y*spins{num}.U(:,:,p-1)*e0) ;     
                    spins{num}.gh(p,3) = 1/opt.epsilon*imag(ez_De*spins{num}.U(:,:,p-1)*e0) ; 

                    spins{num}.gf(p,1) = 1/opt.epsilon*imag(ez_b1x*spins{num}.U(:,:,p-1)*ez) ; 
                    spins{num}.gf(p,2) = 1/opt.epsilon*imag(ez_b1y*spins{num}.U(:,:,p-1)*ez) ;     
                    spins{num}.gf(p,3) = 1/opt.epsilon*imag(ez_De*spins{num}.U(:,:,p-1)*ez) ; 
                else
                    spins{num}.gh(p,1) = 1/opt.epsilon*imag(ez_b1x*e0) ; 
                    spins{num}.gh(p,2) = 1/opt.epsilon*imag(ez_b1y*e0) ;     
                    spins{num}.gh(p,3) = 1/opt.epsilon*imag(ez_De*e0) ; 

                    spins{num}.gf(p,1) = 1/opt.epsilon*imag(ez_b1x*ez) ; 
                    spins{num}.gf(p,2) = 1/opt.epsilon*imag(ez_b1y*ez) ;     
                    spins{num}.gf(p,3) = 1/opt.epsilon*imag(ez_De*ez) ; 
                end


            elseif strcmp(opt.mode,'exact')% B1Phi

                optParam_b1   = W(p,:); 
                del = getDelais_changementvar(opt,W);
                optParam_b1(3) = del(p);
                optParam_b1(1) = optParam_b1(1)*opt.mu;

                [A_1, A_2] = Diff_matrice_de_passage2(spins{num},opt,optParam_b1);
                [A_3] = Diff_matrice_Relax2(spins{num},opt,optParam_b1(3));
                ez_b1x    = ez_it'*getRelaxMat(spins{num},optParam_b1(3),EXPMAT ,opt.mode)*A_1 ;
                ez_b1y    = ez_it'*getRelaxMat(spins{num},optParam_b1(3),EXPMAT ,opt.mode)*A_2 ;
                ez_De     = ez_it'*A_3* getExcMat(opt, spins{num},optParam_b1(1:2)) ;

                ez_it   = (ez_it'*getRelaxMat(spins{num},optParam_b1(3),EXPMAT ,opt.mode) * getExcMat(opt, spins{num},optParam_b1(1:2)))';

                if p>1
                    spins{num}.gh(p,1) = ez_b1x*spins{num}.U(:,:,p-1)*e0 ; 
                    spins{num}.gh(p,2) = (ez_b1y*spins{num}.U(:,:,p-1)*e0) ;     
                    spins{num}.gh(p,3) = (ez_De*spins{num}.U(:,:,p-1)*e0) ; 

                    spins{num}.gf(p,1) = (ez_b1x*spins{num}.U(:,:,p-1)*ez) ; 
                    spins{num}.gf(p,2) = (ez_b1y*spins{num}.U(:,:,p-1)*ez) ;     
                    spins{num}.gf(p,3) = (ez_De*spins{num}.U(:,:,p-1)*ez) ; 
                else
                    spins{num}.gh(p,1) = (ez_b1x*e0) ; 
                    spins{num}.gh(p,2) = (ez_b1y*e0) ;     
                    spins{num}.gh(p,3) = (ez_De*e0) ; 

                    spins{num}.gf(p,1) = (ez_b1x*ez) ; 
                    spins{num}.gf(p,2) = (ez_b1y*ez) ;     
                    spins{num}.gf(p,3) = (ez_De*ez) ; 
                end

            end 
    %      if   strcmp('@getCO',opt.costFunction) 
    %          
    %          temp = 0;
    %          gradS = ((A + lambda*H_U)/((1-lambda*F_U)^2)*lambda*spins{num}.gf(p,:) + lambda/(1-lambda*F_U)*spins{num}.gh(p,:)); % gradient du steady-state 
    %          for j =1:numel(opt.vec)
    %         
    %              if (spins{num}.max == -1)&& (spins{num}.saturation == false)
    % 
    %                      temp = temp + opt.sat*2*signal(j)*gradS*fact(j);
    % 
    % 
    %              elseif  (spins{num}.saturation == false)
    %   
    %                     temp = temp + 2*signal(j)*gradS*fact(j);
    % 
    %              end 
    %          end
    %          spins{num}.g(p,:) = temp;
    %          
              if   strcmp('@getWG_epsilon',opt.costFunction) 

                 temp = 0;
                 gradS = ((A + lambda*H_U)/((1-lambda*F_U)^2)*lambda*spins{num}.gf(p,:) + lambda/(1-lambda*F_U)*spins{num}.gh(p,:)); % gradient du steady-state 
                 D_norme = signal*power(opt.epsilon_abs^2 + signal.^2,-0.5);

                 for j =1:numel(opt.vec)

                     if (spins{num}.max == -1)&& (spins{num}.saturation == false)

                             temp = temp + opt.sat*D_norme(j)*gradS*fact(j);


                     elseif  (spins{num}.saturation == false)

                            temp = temp + D_norme(j)*gradS*fact(j);

                     end 

                 end

                 spins{num}.g(p,:) = temp;

              elseif strcmp('@getWG_epsilon_moins',opt.costFunction)
                 
                 temp = 0;
                 gradS = ((A + lambda*H_U)/((1-lambda*F_U)^2)*lambda*spins{num}.gf(p,:) + lambda/(1-lambda*F_U)*spins{num}.gh(p,:)); % gradient du steady-state 
                 D_norme = signal*power(opt.epsilon_abs^2 + signal.^2,-0.5);

                 for j =1:numel(opt.vec)

                     if (spins{num}.max == -1)&& (spins{num}.saturation == false)

                             temp = temp + opt.sat*D_norme(j)*gradS*fact(j);


                     elseif  (spins{num}.saturation == false)

                            temp = temp + gradS*fact(j);

                     end 

                 end

                 spins{num}.g(p,:) = temp;
                  
              elseif strcmp('@getWG_epsilon_B0',opt.costFunction)
                  temp = 0;
                  gradS = ((A + lambda*H_U)/((1-lambda*F_U)^2)*lambda*spins{num}.gf(p,:) + lambda/(1-lambda*F_U)*spins{num}.gh(p,:)); % gradient du steady-state 
                  for j =1:numel(opt.vec)

                     if (spins{num}.max == -1)&& (spins{num}.saturation == false)
                             temp = temp + opt.sat*gradS*fact(j);
                     elseif  (spins{num}.saturation == false)
                            temp = temp + gradS*fact(j);
                     end 

                  end % boucle sur les
                  spins{num}.g(p,:) = temp;
                  
              end % fonctions de cout
              
        end % boucle sur les impulsions 


    
    end % w0
end % tissues

g= 0;% Compute final gradient

if strcmp('@getWG_epsilon_B0',opt.costFunction)

    for tissue = 1:opt.N_iso % differents tissues
        
        temp = 0;
        for W0 = 1:numel(opt.offsetVecHz)
            
            num = (tissue - 1)*numel(opt.offsetVecHz) + W0;
            if (spins{num}.saturation == false)
                temp =  temp +  spins{num}.g;
            end
            
        end
        grad = somme(tissue)*power(opt.epsilon_abs^2 + somme(tissue).^2,-0.5)*temp;  % somme/sqrt(somme +e2)*Dsomme
        if (spins{num}.max == -1)&& (spins{num}.saturation == false)
            g = g + opt.sat*grad;
        elseif (spins{num}.saturation == false)
            g = g - grad;
        end
    end
    
else

    for num =1:numel(opt.offsetVecHz):numel(spins)   
        if (spins{num}.max == -1)&& (spins{num}.saturation == false)
            g =  g +  spins{num}.g;
        elseif (spins{num}.saturation == false)
            g =  g -  spins{num}.g;
        end
    end
    
end

g(:,end) =  D_ti_alphai'*g(:,end);
 
if strcmp(opt.mode,'exp') % B1xB1y
     
    g(1:opt.Np,1:2) = opt.mu*g(1:opt.Np,1:2);
    
elseif strcmp(opt.mode,'exact')% B1Phi
    
    g(1:opt.Np,1) = opt.mu*g(1:opt.Np,1);
    
end

end 




