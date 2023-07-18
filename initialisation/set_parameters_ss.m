function [opt,samples] = set_parameters_ss(k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Relaxation times of the target tissues in the same order
% T1 : 1450 923 4200 % GM, WM, CSF
% T2 : 96 70 2000
% PD : 0.75 0.65 1

samples.T1 = [1450 923]*10^-3; % GM, WM
samples.T2 = [96 70]*10^-3;
samples.PD = [0.75 0.65];
samples.maxi = [0 1];
samples.saturation = [1 0];

% samples.T1 = [1450 923 4200]*10^-3; % GM, WM, CSF
% samples.T2 = [96 70 2000]*10^-3;
% samples.PD = [0.75 0.65 1]; % Proton density of the target tissues in the same order
% % tissue to be maximized (1) and minimized (0)
% samples.maxi = [1 0 0];
% % tissue to be saturated by constraint (1), (0) otherwise
% samples.saturation = [0 1 1];



% Cost
opt.costFunction = '@getWG_epsilon';%'@getWG_epsilon_moins'
opt.vec = [1]; % Possibilities to take into account multiple lines in the cost function, here if opt.vec = [1] it is only the first line in the train of acquisitions
opt.line_for_constrained_saturation = 1;

opt.saturation_tolerance = 0.0001;  % replace
opt.Np = k+1; % linked to the number of variables to optimize
opt.tp     = 1e-3; % pulse duration (in ms)
opt.expmFunction        = '@expmdemo1';

opt.sat = 1; %  if no strict constraint can put a weight in replace of forced constraint
opt.epsilon_abs = 10^-5; % for wg_epsilon

% Acquisition Parameters
opt.TR = 6e-3; % (s) TR
opt.alpha = 12*pi/180;% (rad) flip angle GRE acquisitions 
opt.Nlignes = 64; %64 lines filled per cycle 
opt.Ncylces = 4; %number of cycles to fill a slide

% Time constraint
opt.time_of_a_segment = (5000)*10^-3; % (s)

%% For Grape
opt.TolX                = 1e-15 ; % Step Tolerance
opt.TolFun              = 1e-5 ; % Function Tolerance
opt.epsilon     = 1e-10 ;
opt.displayIter = 'iter' ;
opt.optimFunction       = 'fmincon';

% Change of variable
opt.changement_of_variable = true ; 
opt.mu = 1000; % b1 = mu * optparam (linear scaling of the optimized b1 values, optparam )

% matrix treatment
opt.mode= 'exp'; %'exp';'exact

% limit
opt.Wmax   = pi/(opt.tp)/opt.mu ;
%to be suppress in the future
opt.TE = 3.2e-3;
opt.TA = 0e-3 ; % ms
opt.TB =0e-3 ; %ms

opt.initVec  = vecteur_initialisation(opt,samples);
opt.gradFunction = '@getGrad_ss'; 
opt.propaFunction = '@propaFunction_ss';
opt.contrainteFunction = '@getcontrainteMz';

if strcmp(opt.mode,'exact')
    opt.propaFunction = '@propafunction_ss';
    opt.initVec = convert2anglephase(opt.initVec,opt);
end

end


