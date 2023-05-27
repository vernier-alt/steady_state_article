function [initVec] = vecteur_initialisation(opt,samples)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


if opt.Np ==1
    
    initVec = zeros(2,3) ;
    initVec(:,3) =140e-3;
    initVec(end,1) = 10*pi/180;
    initVec(end,3) = 4e-3;
    
elseif opt.Np ==2 % Une impulsion
    
%     nombre = 1; 
%     initVec = zeros(1,3,nombre+1) ;
%     angles = [0:1:nombre]/(nombre)*pi/(2*opt.tp);
%     initVec(1,1,:) =reshape( angles,1,1,nombre+1);
% %     initVec(:,1) = pi/(opt.tp);
%     initVec(:,3,:) = 140e-3;
    

    
    N1 = 5; % 10 angles différents
    tab_angles = linspace(90,180,N1)/(opt.tp*opt.mu)*pi/180 ;

    if opt.optimisation_TR
        N2 = 3; % TR differents
        tab_TR = linspace(opt.TRmin,opt.TRmax,N2);
    else
        N2 = 1;
        tab_TR =  [opt.TR];
    end
        
    N3 = 50; % deux emplacements différents pour l'impulsion
    
    if opt.optimisation_TR || opt.optimisation_alpha
        initVec = zeros(3,3,N1*N2*N3) ;
    else
        initVec = zeros(2,3,N1*N2*N3) ;
    end
    
   indice = 0;
    
        for n_angles = 1:N1

            angle = tab_angles(n_angles);

            for n_TR = 1:N2

                TR = tab_TR(n_TR);
                temps_restant = opt.tempsfixe_valeur - opt.Nlignes*TR;
                tab_delais1 = linspace(10e-3,temps_restant-10e-3,N3); % ne pas mettre l'impulsion juste avant l'acquisition

                for n_delais = 1:N3

                    Delais1 = tab_delais1(n_delais);

                    indice = indice +1;

                    % delais
                    initVec(1,3,indice) = Delais1;
                    initVec(2,3,indice) = temps_restant-Delais1;

                    % angles
                    initVec(2,1,indice) = angle;

                    % TR
                    if opt.optimisation_TR
                        initVec(3,3,indice) = TR;
                    end

                    % angle de bascule 
                    if opt.optimisation_alpha
                        initVec(3,1,indice) = min(max(opt.alphamin,1.5*acos(exp(-TR/samples.T1(1)))),opt.alphamax);
                    end
                end 
            end
        end
      size(initVec)
    
elseif opt.Np ==3 
 
    


    
%     N1 = 5; % 5 angles différents pour les impulsions : 0°, 45°, 90°, 135°, 180°
%     N2 = 3; % TR différents
%     N3 = 50; %50% partage du temsp en parts égales (eventuellement p = nchoosek(N3-1,2)) ;
%     p = nchoosek(N3-2,2); 
%     
%     tab_angles = linspace(0,180,N1)/(opt.tp*opt.mu)*pi/180 ;
% 
%     if opt.optimisation_TR
%         tab_TR = linspace(opt.TRmin,opt.TRmax,N2);
%     else
%         N2 = 1;
%         tab_TR =  [opt.TR];
%     end
%     
%     if opt.optimisation_TR || opt.optimisation_alpha
%         initVec = zeros(4,3,N1*N1*N2*p) ;
%     else
%         initVec = zeros(3,3,N1*N1*N2*p) ;
%     end
%     
%     indice = 1 ; 
%     
%     for n_angles1 = 1:N1  % angle 
%         
%         angle1 = tab_angles(n_angles1);
%         
%         for n_angles2 = 1:N1 % angle 
%             
%             angle2 = tab_angles(n_angles2);
%             
%                 for n_TR = 1:N2 % TR 
% 
%                     TR = tab_TR(n_TR);
%                     temps_restant = opt.tempsfixe_valeur - opt.Nlignes*TR;
% 
%                     tab_delais = linspace(1e-3,temps_restant,N3); % ne pas mettre l'impulsion juste avant l'acquisition
% 
% 
%                     for n_delais1 = 2:N3-2 % delais
% 
%                         Delais1 = tab_delais(n_delais1);
% 
%                             for n_delais2 = n_delais1+1:N3-1
% 
%                                 Delais2 = tab_delais(n_delais2);
%                                 
%                                     % delais
%                                     initVec(1,3,indice) = Delais1;
%                                     initVec(2,3,indice) = Delais2-Delais1;
%                                     initVec(3,3,indice) = temps_restant-Delais2;
% 
%                                     if min(initVec(:,3,indice))<0
%                                         break;
%                                     end
% 
%                                     % angles
%                                     initVec(2,1,indice) = angle1;
%                                     initVec(3,1,indice) = angle2;
%                                 
%                                     % TR
%                                     if opt.optimisation_TR
%                                         initVec(opt.Np+1,3,indice) = TR;
%                                     end
% 
%                                     % angle de bascule 
%                                     if opt.optimisation_alpha
%                                         initVec(opt.Np+1,1,indice) = min(max(opt.alphamin,acos(1.001*opt.Kmin*exp(+TR/samples.T1(1)))),opt.alphamax);
%                                     end
% 
%                                     indice = indice + 1; 
%                                 
%                             end
% 
%                     end
%                  end 
%         end
%     end
%      
% 
%   indice
%   size(initVec  )
%   
%   disp('stop')
  
  %% 2 Impulsions 
  
    N1 = 50;
    N2 = 3; % TR différents
%     if opt.optimisation_TR
%         tab_TR = linspace(opt.TRmin,opt.TRmax,N2);
%     else
        N2 = 1;
        tab_TR =  [opt.TR];
%     end

%     if opt.optimisation_TR || opt.optimisation_alpha
%         initVec = zeros(4,3,100) ;
%     else
        initVec = zeros(3,3,100) ;
%     end
  
  indice = 1;
  
  
  for num_Tr = 1:numel(tab_TR)
  
      TR = tab_TR(num_Tr);
      temps_restant = opt.tempsfixe_valeur - opt.Nlignes*TR;
      
      % 0 impulsion
      % delais
        initVec(1,3,indice) = temps_restant/3;
        initVec(2,3,indice) = temps_restant/3;
        initVec(3,3,indice) = temps_restant/3;

        % angles
        initVec(2,1,indice) = 0;
        initVec(3,1,indice) = 0;

        % TR
%         if opt.optimisation_TR
%             initVec(opt.Np+1,3,indice) = TR;
%         end

        % angle de bascule 
%         if opt.optimisation_alpha
%             initVec(opt.Np+1,1,indice) = min(max(opt.alphamin,acos(1.001*opt.Kmin*exp(+TR/samples.T1(1)))),opt.alphamax);
%         end
        
        indice = indice + 1;

      % 1 impulsion : 180 ou 90 à N1 emplacements
        temps = linspace(0,temps_restant,N1+2);
        w1x = [90 180]/(opt.tp*opt.mu)*pi/180;
        w1y = [0 0]/(opt.tp*opt.mu)*pi/180;
        for w = 1:2
            for tps = 2:N1
                initVec(1,3,indice) = temps(tps);
                initVec(2,3,indice) = (temps_restant-temps(tps))/2 ;
                initVec(3,3,indice) = (temps_restant-temps(tps))/2;

                % angles
                initVec(2,1,indice) = w1x(w);
                initVec(2,2,indice) = w1y(w);

                % TR
%                 if opt.optimisation_TR
%                     initVec(opt.Np+1,3,indice) = TR;
%                 end

%                 % angle de bascule 
%                 if opt.optimisation_alpha
%                     initVec(opt.Np+1,1,indice) = min(max(opt.alphamin,acos(1.001*opt.Kmin*exp(+TR/samples.T1(1)))),opt.alphamax);
%                 end
                indice = indice + 1;
            end
        end

        % 2 impulsions 

        % si 1ere impulsion = 90° (-90° serait pareil) N1 possibilités pour la placer, suivi avec un délais de
        % max 20ms par soit une 90°, soit une -90° , soit 180°
        % à N1 emplacements différents

        temps = linspace(0,temps_restant,N1+2);
        w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
        w2y = [0 0 0]/(opt.tp*opt.mu)*pi/180;

        for w = 1:3
            for tps = 2:N1
                initVec(1,3,indice) = temps(tps);
                initVec(2,3,indice) = min(20e-3,temps_restant/(N1+2)) ;
                initVec(3,3,indice) = temps_restant-initVec(2,3,indice);

                % angles
                initVec(2,1,indice) = 90/(opt.tp*opt.mu)*pi/180;% 90°
                initVec(2,2,indice) = 0;
                initVec(3,1,indice) = w2x(w); %90 ou -90° ou 180°
                initVec(3,2,indice) = w2y(w);

                % TR
%                 if opt.optimisation_TR
%                     initVec(opt.Np+1,3,indice) = TR;
%                 end

                % angle de bascule 
%                 if opt.optimisation_alpha
%                     initVec(opt.Np+1,1,indice) = min(max(opt.alphamin,acos(1.001*opt.Kmin*exp(+TR/samples.T1(1)))),opt.alphamax);
%                 end
                indice = indice + 1;
            end
        end

        % si 1ere impulsion = 180°, N1 possibilités pour la placer
        % la seconde peut être une 90° ou -90° ou 180° sur les emplacements
        % restants
        temps = linspace(0,temps_restant,N1+2);
        w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
        w2y = [0 0 0]/(opt.tp*opt.mu)*pi/180;

        for w = 1:3
            for tps = 2:N1
                for tps2 = tps+1:N1

                    initVec(1,3,indice) = temps(tps);
                    initVec(2,3,indice) = temps(tps2)-temps(tps) ;
                    initVec(3,3,indice) = temps_restant-initVec(2,3,indice);

                    % angles
                    initVec(2,1,indice) = 180/(opt.tp*opt.mu)*pi/180; % 180°
                    initVec(2,2,indice) = 0;
                    initVec(3,1,indice) = w2x(w); % 180°, 90°, -90°
                    initVec(3,2,indice) = w2y(w);

                    % TR
%                     if opt.optimisation_TR
%                         initVec(opt.Np+1,3,indice) = TR;
%                     end

                    % angle de bascule 
%                     if opt.optimisation_alpha
%                         initVec(opt.Np+1,1,indice) = min(max(opt.alphamin,acos(1.001*opt.Kmin*exp(+TR/samples.T1(1)))),opt.alphamax);
%                     end
                    indice = indice + 1;
                end
            end
        end
  end   
indice



elseif opt.Np == 4

  %% 3 Impulsions 
  
    N1 = 20;
    initVec = zeros(4,3,100) ;
    indice = 1;
    temps_restant = opt.tempsfixe_valeur - opt.Nlignes*opt.TR;
    % 0 impulsion
    % delais
    initVec(1,3,indice) = temps_restant/4;
    initVec(2,3,indice) = temps_restant/4;
    initVec(3,3,indice) = temps_restant/4;
    initVec(4,3,indice) = temps_restant/4;
    % angles
    initVec(2,1,indice) = 0;
    initVec(3,1,indice) = 0;
    initVec(4,1,indice) = 0;
    % 1 impulsion : 180, 45, 90 à N1 emplacements
    temps = linspace(0,temps_restant,N1+2);
    w1x = [90 45 180]/(opt.tp*opt.mu)*pi/180;
    w1y = [0 0 0]/(opt.tp*opt.mu)*pi/180;
        for w = 1:length(w1x)
            for tps = 2:N1
                
                % delays
                initVec(1,3,indice) = temps(tps);
                initVec(2,3,indice) = (temps_restant-temps(tps))/3 ;
                initVec(3,3,indice) = (temps_restant-temps(tps))/3 ;
                initVec(4,3,indice) = (temps_restant-temps(tps))/3 ;

                % angles
                initVec(2,1,indice) = w1x(w);
                initVec(2,2,indice) = w1y(w);
                indice = indice + 1;
            end
        end

        % 2 impulsions 
        % si 1ere impulsion = 90° (-90° serait pareil) N1 possibilités pour la placer, suivi avec un délais de
        % max 20ms par soit une 90°, soit une -90° , soit 180°
        % à N1 emplacements différents

        temps = linspace(0,temps_restant,N1+2);
        w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
        w2y = [0 0 0 ]/(opt.tp*opt.mu)*pi/180;

        for w = 1:length(w2x)
            for tps = 2:N1
                initVec(1,3,indice) = temps(tps);
                initVec(2,3,indice) = min(20e-3,temps_restant/(N1+2)) ;
                initVec(3,3,indice) = (temps_restant-initVec(2,3,indice))/2;
                initVec(4,3,indice) = (temps_restant-initVec(2,3,indice))/2;

                % angles
                initVec(2,1,indice) = 90/(opt.tp*opt.mu)*pi/180;% 90°
                initVec(2,2,indice) = 0;
                initVec(3,1,indice) = w2x(w); %90 ou -90° ou 180°
                initVec(3,2,indice) = w2y(w);
                
                indice = indice + 1;
            end
        end

        % si 1ere impulsion = 180°, N1 possibilités pour la placer
        % la seconde peut êopt.TRe une 90° ou -90° ou 180° sur les emplacements
        % restants
        temps = linspace(0,temps_restant,N1+2);
        w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
        w2y = [0 0 0]/(opt.tp*opt.mu)*pi/180;

        for w = 1:3
            for tps = 2:N1-1
                for tps2 = tps+1:N1

                    initVec(1,3,indice) = temps(tps);
                    initVec(2,3,indice) = temps(tps2)-temps(tps) ;
                    initVec(3,3,indice) = (temps_restant-initVec(2,3,indice))/2;
                    initVec(4,3,indice) = (temps_restant-initVec(2,3,indice))/2;

                    % angles
                    initVec(2,1,indice) = 180/(opt.tp*opt.mu)*pi/180; % 180°
                    initVec(2,2,indice) = 0;
                    initVec(3,1,indice) = w2x(w); % 180°, 90°, -90°
                    initVec(3,2,indice) = w2y(w);

                    indice = indice + 1;
                end
            end
        end
        
        % 3 impulsions
        % on reprend le schéma précédent
        % si premiere impulsions 90 - 180 , la opt.troisieme suit avec un délais
        % de 20ms (car toujours aimantation opt.TRansverse)
        % choix parmi 90 ou -90
        temps = linspace(0,temps_restant,N1+2);
        w3x = [90 -90 ]/(opt.tp*opt.mu)*pi/180;
        w3y = [0 0]/(opt.tp*opt.mu)*pi/180;

        for w = 1:2
            for tps = 2:N1
                initVec(1,3,indice) = temps(tps);
                initVec(2,3,indice) = min(20e-3,temps_restant/(N1+2)) ;
                initVec(3,3,indice) = min(20e-3,temps_restant/(N1+2));
                initVec(4,3,indice) = (temps_restant-initVec(2,3,indice)-initVec(3,3,indice));


                % angles
                initVec(2,1,indice) = 90/(opt.tp*opt.mu)*pi/180;% 90°
                initVec(2,2,indice) = 0;
                initVec(3,1,indice) = 180/(opt.tp*opt.mu)*pi/180;
                initVec(3,2,indice) = 0;
                initVec(4,1,indice) = w3x(w); %90 ou -90° 
                initVec(4,2,indice) = w3y(w);

                indice = indice + 1;
            end
        end
        
        % si premiere impulsions 90 - -90  / -90 - 90 ce qui est la même
        % chose ou 90 - 90
        % est une 180, placé sur les N1-2 places restantes 
        temps = linspace(0,temps_restant,N1+2);

        w2x = [90 -90 ]/(opt.tp*opt.mu)*pi/180;
        w2y = [0 0 ]/(opt.tp*opt.mu)*pi/180;
        
        for w = 1:2
            for tps = 2:N1-1
                for tps2 = tps+1:N1
                    initVec(1,3,indice) = temps(tps);
                    initVec(2,3,indice) = min(20e-3,temps_restant/(N1+2)) ;
                    initVec(3,3,indice) = temps(tps2)-initVec(2,3,indice)-initVec(1,3,indice);
                    initVec(4,3,indice) = temps_restant-sum(initVec(1:3,3,indice));


                    % angles
                    initVec(2,1,indice) = 90/(opt.tp*opt.mu)*pi/180;% 90°
                    initVec(2,2,indice) = 0;
                    initVec(3,1,indice) = w2x(w);%90 ou -90
                    initVec(3,2,indice) = w2y(w);
                    initVec(4,1,indice) = 180/(opt.tp*opt.mu)*pi/180;% 180°
                    initVec(4,2,indice) = 0;

                    indice = indice + 1;
                end
            end
        end
        
        % si premiere impulsions 180 - 180
        % choix : 180
        % est une 180, placé sur les N1-3 places restantes 
        temps = linspace(0,temps_restant,N1+2);

        for w = 1:3
            for tps = 2:N1-2
                for tps2 = tps+1:N1-1
                    for tps3 = tps2+1:N1
                        initVec(1,3,indice) = temps(tps);
                        initVec(2,3,indice) = temps(tps2)-temps(tps) ;
                        initVec(3,3,indice) = temps(tps3)-sum(initVec(1:2,3,indice));
                        initVec(4,3,indice) = temps_restant-sum(initVec(1:3,3,indice));


                        % angles
                        initVec(2,1,indice) = 180/(opt.tp*opt.mu)*pi/180;% 180°
                        initVec(2,2,indice) = 0;
                        initVec(3,1,indice) = 180/(opt.tp*opt.mu)*pi/180;% 180°
                        initVec(3,2,indice) = 0;
                        initVec(4,1,indice) = 180/(opt.tp*opt.mu)*pi/180;% 180°
                        initVec(4,2,indice) = 0;

                        indice = indice + 1;
                    end
                end
            end
        end
  end      
    
    
indice
    
end







