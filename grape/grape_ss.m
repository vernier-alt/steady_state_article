function [historyOpt,spins]  = grape_ss(spins,opt) 
% W is the optimal control field, fval the correspding cost
% % Set up shared variables with OUTFUN
historyOpt.W    = [];
historyOpt.fval = [];
DYN             = str2func( opt.propaFunction ) ;
GRAD = str2func( opt.gradFunction ) ;
  
% format long 
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
        'GradObj','on',... // specify objective gradient 
        'Gradconstr','on',...// specify constraint gradient 
        'Hessian',{'lbfgs',5},...  
        'DerivativeCheck','on',...
        'AlwaysHonorConstraints','none');  
    
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
        
    % No pulse at the beginning of the segment/cycle 
    lbound(1,1:2) = 0;
    ubound(1,1:2) = 0;

    % Non linear constraint to force to null the signal of a tissue:
    if opt.contrainte
        nonlcon = @contrainte;
    else
        nonlcon = [];
    end
    
    % Linear constraints to respect that the sum of the delays is equal to
    % the segment duration :
    % opt.W_ini is a mxn matrix, it is reorganized by the solver as a
    % vector (let's call it W)
    % of mxn lines, each column is put under the previous one
    % linear constraint is : A*W = 0 or A*W <=0 if tolerance
    Tmin = 70e-3; % minimum duration for the last delay of the preparation
    if opt.changement_of_variable
        ll = size(opt.W_ini,1);
        A = zeros(1,ll*3);
        b = zeros(1,1);
        A(1,end) = -(opt.time_of_a_segment-opt.Nlignes*opt.TR) + Tmin;
        A(1,ll*2+1:end-1) = Tmin;
        lbound(1:opt.Np,3) = 1e-3; % Optimized variables must be positive, to be too much close to zero may cause numerical problems (divisions in the conversion)
        ubound(1:opt.Np,3) = inf;
        W = fmincon(@getCost,opt.W_ini,A,b,[],[],lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % remain an inequlity on the optimized variables as the last delay of the preparation must be superior to Tmin
    else
        lbound(1:opt.Np,3) = 0; % Optimized delays must be positive
        lbound(opt.Np,3) = Tmin;
        ubound(1:opt.Np,3) = inf;
        A = zeros(1,size(opt.W_ini,1)*3);
        A(1,end-opt.Np+1:end) = 1;
        b = opt.time_of_a_segment-opt.Nlignes*opt.TR ;
        W = fmincon(@getCost,opt.W_ini,[],[],A,b,lbound,ubound,nonlcon,optionsFmincon,spins,opt)   ; % equality on the sum of the delay 
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
    % Compute the transfer matrix of the preparation
    for si = 1:numel(spins)
        [spins{si}.U]   = DYN(W,spins{si},opt) ;  %propa function 
    end
    % compute the objective (cost) function
    COST    = str2func( opt.costFunction );
    f       = COST(W,spins,opt) ;
    % compute the gradient of the objective function
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
      % Compute the constraint function for each spin and their gradients
      % (stored in a matrix)
      % can be inegality of strict equality : 
      c = [];
      Dc = [];
      ceq = [];
      Dceq = [];
      for si = 1:numel(spins)
        [spins{si}.U]   = DYN(W,spins{si},opt) ;  %propa function requires W to be in the temporal domain
      end
      opt.line_for_constrained_saturation = 1 ; % contrainte sur la premiere ligne
      CONT    = str2func( opt.contrainteFunction );    
      f = opt.line_for_constrained_saturation;
      if opt.saturation_tolerance>0
          c = CONT(W,spins,opt,f) - opt.saturation_tolerance;
          Dc = getGrad_constrainte_ss(W,spins,opt,f);
      else
          ceq = CONT(W,spins,opt,f);
          Dceq = getGrad_constrainte_ss(W,spins,opt,f);
      end
  
end

end
