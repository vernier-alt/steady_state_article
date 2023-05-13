function [historyOpt,spins]  = grape_ss(spins,opt) 
% W is the optimal control field, fval the correspding cost
% % Set up shared variables with OUTFUN
historyOpt.W    = [];
historyOpt.fval = [];
DYN             = str2func( opt.propaFunction ) ;
GRAD = str2func( opt.gradFunction ) ;

% format long :   ezezefzerf efezf  vhjvv
if strcmp(opt.optimFunction,'fmincon')    

    optionsFmincon = optimoptions(@fmincon,'OutputFcn',@outfun,...    
        'TolX',opt.TolX,...
        'TolFun',opt.TolFun,...
        'Display',opt.displayIter,...
        'MaxIter',opt.maxIter,...
        'MaxFunEvals',Inf,...
        'Algorithm','sqp',...
        'InitBarrierParam',0.5,...
        'InitTrustRegionRadius',10,...        
        'GradObj','on',...
        'Gradconstr','on',...
        'Hessian',{'lbfgs',5},...  
        'DerivativeCheck','off',...
        'AlwaysHonorConstraints','none');  
    
    %%'Gradconstr','off',...
    
    
    lbound = zeros(size(opt.W_ini)) ; %lower amplitude constraint 
    ubound = zeros(size(opt.W_ini)) ; %upper amplitude constraint 
    
    if strcmp(opt.mode,'exp')
        lbound(1:opt.Np,1:2) = -1.001*opt.Wmax;   % -opt.Wmax ; %lower amplitude constraint 
        ubound(1:opt.Np,1:2)  = 1.001*opt.Wmax ; %upper amplitude constraint 
    else
        lbound(1:opt.Np,1) =  0;    
        ubound(1:opt.Np,1)  = 1.001*opt.Wmax ; 
        lbound(1:opt.Np,2) =  0;    
        ubound(1:opt.Np,2)  = pi ; 
    end    
        
    % Angle de bascule
    if  opt.optimisation_alpha
        lbound(opt.Np+1,1) = opt.alphamin;
        ubound(opt.Np+1,1) = opt.alphamax;
    end

    
    % Première impulsion nulle 
    lbound(1,1:2) = 0;
    ubound(1,1:2) = 0;
    
%     % Pas d'impulsion selon z
%     lbound(1:end,2) = 0;
%     ubound(1:end,2) = 0;
%  
    % Force t2-prep
%     lbound(2,1) = pi/(2*opt.tp)/opt.mu;
%     ubound(2,1) = pi/(2*opt.tp)/opt.mu;
%     lbound(3,1) = pi/(2*opt.tp)/opt.mu;
%     ubound(3,1) = pi/(2*opt.tp)/opt.mu;
%     lbound(4,1) = pi/(opt.tp)/opt.mu;
%     ubound(4,1) = pi/(opt.tp)/opt.mu;
%     lbound(2:4,2) = 0;
%     ubound(2:4,2) = 0;

    % Force DIR
%     lbound(2,1) = pi/(opt.tp)/opt.mu;
%     ubound(2,1) = pi/(opt.tp)/opt.mu;
%     lbound(3,1) = pi/(opt.tp)/opt.mu;
%     ubound(3,1) = pi/(opt.tp)/opt.mu;
%     lbound(2:3,2) = 0;
%     ubound(2:3,2) = 0;



%      W = fmincon(@getCost,opt.W_ini,A,b,[],[],lbound,ubound,[],optionsFmincon,spins,opt)   ;
if opt.contrainte
    nonlcon = @contrainte;
else
    nonlcon = [];
end

if opt.optimisation_TR
    if opt.changement_variable_delais 
        % alpha - TR
        lbound(1:end,3) = 1e-3;
        ubound(1:end,3) = inf;
        A = zeros(2,(opt.Np+1)*3);
        b = zeros(2,1);
        A(2,end) = -opt.tempsfixe_valeur + opt.Nlignes*opt.TRmin;
        A(2,end-3:end-1) = opt.TRmin;
        A(1,end) = opt.tempsfixe_valeur - opt.Nlignes*opt.TRmax;
        A(1,end-3:end-1) = -opt.TRmax;
        disp('contrainte')
        W = fmincon(@getCost,opt.W_ini,A,b,[],[],lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % inigalité pour respecter le min TR
    else
        lbound(end,3) = opt.TRmin;
        ubound(end,3) = opt.TRmax;
        A = zeros(1,(opt.Np+1)*3);
        A(1,end-(opt.Np+1)+1:end-1) = 1;
        A(1,end) = opt.Nlignes;
        b = [opt.tempsfixe_valeur] ;
        W = fmincon(@getCost,opt.W_ini,[],[],A,b,lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % égalite et borne min max TR 
    end
else
    if opt.changement_variable_delais 
        ll = size(opt.W_ini,1);
        A = zeros(2,ll*3);
        b = zeros(2,1);
        Tmin = 10e-3;
        A(1,end) = -(opt.tempsfixe_valeur-opt.Nlignes*opt.TR) + Tmin;
        A(1,ll*2+1:end-1) = Tmin;
        
%         A(2,end-ll+1:end) = 10e-3;
%         A(2,end-ll+2) = -(opt.tempsfixe_valeur-opt.Nlignes*opt.TR) + 10e-3;
        
%         Delais_min = 10e-3;
%         Delais_min_dernier = 50e-3;
%         
%         amin = [ones(opt.Np-1,1)*Delais_min; Delais_min_dernier];
%         v = diag(amin')*ones(opt.Np);
%         v = v - eye(opt.Np)*(opt.tempsfixe_valeur-opt.Nlignes*opt.TR);
%         
%         A = zeros(opt.Np,ll*3);
%         b = zeros(opt.Np,1);
%         A(:,ll*2+1:ll*2+opt.Np)= v;

        
        lbound(1:opt.Np,3) = 1e-3;
        ubound(1:opt.Np,3) = inf;
        
        W = fmincon(@getCost,opt.W_ini,A,b,[],[],lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % inigalité pour respecter le min TR
    else
        lbound(1:opt.Np,3) = 0;
        lbound(opt.Np,3) = 50e-3;
        ubound(1:opt.Np,3) = inf;
        
        A = zeros(1,size(opt.W_ini,1)*3);
        A(1,end-opt.Np+1:end) = 1;
        b = opt.tempsfixe_valeur-opt.Nlignes*opt.TR ;
        
        W = fmincon(@getCost,opt.W_ini,[],[],A,b,lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % égalite et borne min max TR 
    end
end
    
end    


% final_Wt = getTemporalControl(W,opt) ;
historyOpt.W = W;
for s = 1:numel(spins)
    [spins{s}.U] = DYN(W,spins{s},opt) ;
end

%% Additional functions
% %----------------------------

function [f,g] = getCost(W,spins,opt)
    % Propagate all spins from M(t=0) to M(t=tf) with the selected propagation
    % function
    
    for si = 1:numel(spins)
        [spins{si}.U]   = DYN(W,spins{si},opt) ;  %propa function requires W to be in the temporal domain
    end
    
    COST    = str2func( opt.costFunction );
    
    f       = COST(W,spins,opt) ;
    g = GRAD(W,spins,opt);
      
end

function stop = outfun(x,optimValues,state,spins,opt)
     stop = false;
 
     switch state
         case 'init'
%              hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           historyOpt.fval  = [historyOpt.fval; optimValues.fval];
           historyOpt.W     = [historyOpt.W; x];%%constrviolation
           historyOpt.constraint = [optimValues.constrviolation];
         case 'done'
%              hold off 
         otherwise
     end
end

function [c , ceq, Dc, Dceq] = contrainte(W,spins,opt) 
    
      % default 
      c = [];
      Dc = [];
      ceq = [];
      Dceq = [];
    
    
      for si = 1:numel(spins)
        [spins{si}.U]   = DYN(W,spins{si},opt) ;  %propa function requires W to be in the temporal domain
      end
    
      f = 1 ; % contrainte sur la premiere ligne
      CONT    = str2func( opt.contrainteFunction );    

      if opt.saturation_tolerance>0
          c = CONT(W,spins,opt,f) - opt.saturation_tolerance;
          Dc = getGrad_constrainte_ss(W,spins,opt,f);
      else
          ceq = CONT(W,spins,opt,f);
          Dceq = getGrad_constrainte_ss(W,spins,opt,f);
      end

      if (opt.optimisation_TR == true) || (opt.optimisation_alpha == true)
          l = size(c,2);
          c = horzcat(c,contrainteT1etoile(W,spins,opt));
          Dc = horzcat(Dc,getGrad_contrainteT1etoile(W,spins,opt));
      end


      
end

end
