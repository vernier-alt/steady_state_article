function Jcontrainte = getGrad_constrainte_ss(W,spins,opt,num,f)

format long 

EXPMAT  = str2func(opt.expmFunction) ;   

ez = [0 0 0 1]';
e0 = [1 0 0 0]';

% sum = 0;

if opt.changement_of_variable

    v = [W(1:end,3)'];
    M1 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
    M1(1:opt.Np,1:opt.Np) = eye(opt.Np)*sum(v);
    M2 = zeros(size(W,1));% au cas ou on optimiserait alpha et pas TR
    M2(1:opt.Np,1:opt.Np) = diag(W(1:opt.Np,3)')*ones(opt.Np);
    D_ti_alphai = (M1-M2)/(sum(v)^2)*(opt.time_of_a_segment-opt.TR*opt.Nlignes);

else
    D_ti_alphai = eye(size(W,1));
end


alpha = opt.alpha;
TR = opt.TR;

for num =1:numel(spins)

 
     if spins{num}.saturation == true
        
        ez_it = [0 0 0 1]';


        EA = exp(-opt.TA/spins{num}.T1); %ms
        EB = exp(-opt.TB/spins{num}.T1); %ms


        E1 = exp(-TR/spins{num}.T1); %ms

        K = cos(alpha)*E1; %ms .

        A = spins{num}.Mt0(end)*((1-EB)+( (1-E1)*(1-K^opt.Nlignes)/(1-K)+K^opt.Nlignes*(1-EA))*EB);
        lambda = K^opt.Nlignes*EA*EB;

        H_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*e0;
        F_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*[0 0 0 1]';

        s = spins{num}.Mt0(end)*(1-E1)/(1-K);
        
        ss = ((A + lambda*H_U)/(1-lambda*F_U ));

        ss_EB = (ss -(1-EB))/EB;
        
        signal = (s +  ( ss_EB -s)*K^(-opt.Nlignes));
        
        fact = (EB^-1)*K^(-opt.Nlignes);

        % dommage de refaire tout �a, a voir 
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

         if strcmp('@getcontrainteMz',opt.contrainteFunction)

                 gradS = ((A + lambda*H_U)/((1-lambda*F_U)^2)*lambda*spins{num}.gf(p,:) + lambda/(1-lambda*F_U)*spins{num}.gh(p,:)); % gradient du steady-state 

                 spins{num}.gc(p,:) = 2*signal*gradS*fact;

         end
         

         end

     end
        
end 



% Jcontrainte = 0 ; % Jacobienne (possiblement plusieurs contraintes )
j = 1;
for num = 1:numel(spins)
    if spins{num}.saturation == true 
         spins{num}.gc(:,end) =  D_ti_alphai'*spins{num}.gc(:,end);
         if strcmp(opt.mode,'exp') % B1xB1y
            spins{num}.gc(1:opt.Np,1:2) = opt.mu*spins{num}.gc(1:opt.Np,1:2);
         elseif strcmp(opt.mode,'exact')% B1Phi
            spins{num}.gc(1:opt.Np,1) = opt.mu*spins{num}.gc(1:opt.Np,1);
         end
         Jcontrainte(:,j) = spins{num}.gc(:);
         j=j+1;
    end
end



%%


end 
