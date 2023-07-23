function [initVec] = vecteur_initialisation(opt,samples)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


if opt.Np ==1
    
    initVec = zeros(2,3) ;
    initVec(:,3) =140e-3;
    initVec(end,1) = 10*pi/180;
    initVec(end,3) = 4e-3;
    
elseif opt.Np ==2 % Une impulsion
    
    N1 = 50; % 10 differents angles
    tab_angles = linspace(90,180,N1)/(opt.tp*opt.mu)*pi/180 ;
    N3 = 50; % 50 different places
    initVec = zeros(2,3,N1*N3) ;
    indice = 1;
    
        for n_angles = 1:N1

            angle = tab_angles(n_angles);

                TR = opt.TR;
                temps_restant = opt.time_of_a_segment - opt.Nlignes*TR;
                tab_delais1 = linspace(10e-3,temps_restant-10e-3,N3); % ne pas mettre l'impulsion juste avant l'acquisition

                for n_delais = 1:N3

                    Delay = tab_delais1(n_delais);

                    

                    % delais
                    initVec(1,3,indice) = Delay;
                    initVec(2,3,indice) = temps_restant-Delay;

                    % angles
                    initVec(2,1,indice) = angle;
                    
                    indice = indice +1;

                end 

        end
      size(initVec)
    
elseif opt.Np ==3 
  % 2 Impulsions 
  
   N1 = 100;

   N2 = 1;
   tab_TR =  [opt.TR];

    initVec = zeros(3,3,100) ;

  
  indice = 1;
  
  
  for num_Tr = 1:numel(tab_TR)
  
      TR = tab_TR(num_Tr);
      temps_restant = opt.time_of_a_segment - opt.Nlignes*TR;
      
      % 0 pulse :
      % delay
        initVec(1,3,indice) = temps_restant/3;
        initVec(2,3,indice) = temps_restant/3;
        initVec(3,3,indice) = temps_restant/3;

        % angles
        initVec(2,1,indice) = 0;
        initVec(3,1,indice) = 0;

        indice = indice + 1;

      % 1 pulse : 45, 180 or 90 at N1 places
        temps = linspace(0,temps_restant,N1+2);
        w1x = [45 90 180]/(opt.tp*opt.mu)*pi/180;
        w1y = [0]/(opt.tp*opt.mu)*pi/180;
        for wy = 1:numel(w1y)
            for wx = 1:numel(w1x)
                for tps1 = 1:N1
                        
                        % delay
                        initVec(1,3,indice) = temps(tps1);
                        initVec(2,3,indice) = (temps_restant-temps(tps1))/2;
                        initVec(3,3,indice) = (temps_restant-temps(tps1))/2;

                        % angles
                        initVec(2,1,indice) = w1x(wx);
                        initVec(2,2,indice) = w1y(wy);

                        indice = indice + 1;
                        

                end
            end
        end

        w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
        w2y = [0]/(opt.tp*opt.mu)*pi/180;
        

        for w1x_it = 1:numel(w1x)
            for w2y_it = 1:numel(w2y)
                for w2x_it = 1:numel(w2x)
                    for tps1 = 2:N1
                        for tps2 = tps1+1:N1+1
                        
                        % delay
                        initVec(1,3,indice) = temps(tps1);
                        initVec(2,3,indice) = temps(tps2)-temps(tps1);
                        initVec(3,3,indice) = (temps_restant-temps(tps2));

                        % angles
                        initVec(2,1,indice) = w1x(w1x_it);
                        initVec(2,2,indice) = 0;
                        initVec(3,1,indice) = w2x(w2x_it); 
                        initVec(3,2,indice) = w2y(w2y_it);
                        
                        indice = indice + 1;
                        
                        end 
                    end
                end
            end 
        end

  end   

elseif opt.Np == 4

  %% 3 Impulsions 
  
    N1 = 30;
    
    initVec = zeros(4,3,100) ;
    indice = 1;
    temps_restant = opt.time_of_a_segment - opt.Nlignes*opt.TR;
    
    % 0 pulse
    % delay
    initVec(1,3,indice) = temps_restant/4;
    initVec(2,3,indice) = temps_restant/4;
    initVec(3,3,indice) = temps_restant/4;
    initVec(4,3,indice) = temps_restant/4;
    % angles
    initVec(2,1,indice) = 0;
    initVec(3,1,indice) = 0;
    initVec(4,1,indice) = 0;
    indice =indice +1;
    
    % 1 pulse : 45, 180 or 90 at N1 places
    temps = linspace(0,temps_restant,N1+3);
    w1x = [45 90 180]/(opt.tp*opt.mu)*pi/180;
    w1y = [0]/(opt.tp*opt.mu)*pi/180;
    for wy = 1:numel(w1y)
        for wx = 1:numel(w1x)
            for tps1 = 1:N1

                    % delay
                    initVec(1,3,indice) = temps(tps1);
                    initVec(2,3,indice) = (temps_restant-temps(tps1))/3;
                    initVec(3,3,indice) = (temps_restant-temps(tps1))/3;
                    initVec(4,3,indice) = (temps_restant-temps(tps1))/3 ;
                    
                    % angles
                    initVec(2,1,indice) = w1x(wx);
                    initVec(2,2,indice) = w1y(wy);
                    initVec(3:4,1:2,indice) = 0;
                    
                    indice = indice + 1;


            end
        end
    end
    
    % 2 pulses : 
    w2x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
    w2y = [0]/(opt.tp*opt.mu)*pi/180;
    for w1x_it = 1:numel(w1x)
        for w2y_it = 1:numel(w2y)
            for w2x_it = 1:numel(w2x)
                for tps1 = 2:N1
                    for tps2 = tps1+1:N1+1

                    % delay
                    initVec(1,3,indice) = temps(tps1);
                    initVec(2,3,indice) = temps(tps2)-temps(tps1);
                    initVec(3,3,indice) = (temps_restant-temps(tps2))/2;
                    initVec(4,3,indice) = (temps_restant-temps(tps2))/2;

                    % angles
                    initVec(2,1,indice) = w1x(w1x_it);
                    initVec(2,2,indice) = 0;
                    initVec(3,1,indice) = w2x(w2x_it); 
                    initVec(3,2,indice) = w2y(w2y_it);
                    initVec(4,1:2,indice) = 0;

                    indice = indice + 1;

                    end 
                end
            end
        end 
    end

    % 3 pulses : 
    w3x = [90 -90 180]/(opt.tp*opt.mu)*pi/180;
    for w1x_it = 1:numel(w1x)
        for w2y_it = 1:numel(w2y)
            for w2x_it = 1:numel(w2x)
                for w3x_it = 1:numel(w3x)
                    for tps1 = 1:N1
                        for tps2 = tps1+1:N1+1
                            for tps3 = tps2+1:N1+2
                                
                                % delay
                                initVec(1,3,indice) = temps(tps1);
                                initVec(2,3,indice) = temps(tps2)-temps(tps1);
                                initVec(3,3,indice) = temps(tps3)-temps(tps2);
                                initVec(4,3,indice) = temps_restant-temps(tps3);

                                % angles
                                initVec(2,1,indice) = w1x(w1x_it);
                                initVec(2,2,indice) = 0;
                                initVec(3,1,indice) = w2x(w2x_it); 
                                initVec(3,2,indice) = 0;
                                initVec(4,1,indice) = w3x(w3x_it);
                                initVec(4,2,indice) = 0;

                                indice = indice + 1;
                        
                            end

                        end 
                    end
                end
            end
        end 
    end
        
    
    
indice
    
end







