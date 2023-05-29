function [opt,samples] = set_parameters_ss(k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Relaxation times of the target tissues in the same order
samples.T1 = [1500 1800]*10^-3; 
samples.T2 = [27 34]*10^-3;
% Proton density of the target tissues in the same order
samples.PD = [1 1];
% tissue to be maximized (1) and minimized (0)
samples.maxi = [1 0];
% tissue to be saturated by constraint (1), (0) otherwise
samples.saturation = [0 1];
% Cost
opt.costFunction = '@getWG_epsilon_moins';%'@getWG_epsilon''@getWG'; % '@getCO'

opt.saturation_tolerance = 0.0001; % tolerance 0.001 c'est bien
opt.Np = k+1; % linked to the number of variables to opimize
opt.tp     = 1e-3;%1e-3; % pulse duration (in ms)
opt.expmFunction        = '@expmdemo1';

opt.sat = 1; % pour getCO, if no strict constraint can put a weight
opt.epsilon_abs = 10^-5; % pour wg_epsilon

%% For Grape
opt.TolX                = 1e-15 ; % Step Tolerance
opt.TolFun              = 1e-5 ; % Function Tolerance
opt.epsilon     = 1e-10 ;
opt.displayIter = 'iter' ;
opt.optimFunction       = 'fmincon';

% Change of variable
opt.changement_variable_delais = true ; % si changement de variable, on conserve une contrainte d'inégalité si TR min
opt.mu = 1000; % b = mu * optparam, si mu trop grand, optparam trop petit, gradient élevé ereur de derivative check

% Les paramètres d'acquisition
opt.TR = 6e-3; % Sinon voici le TR
opt.alpha = 13*pi/180;% Sinon voici le alpha 
opt.Nlignes = 64; %64 % lines per cycle to fill a slice
opt.Ncylces = 4;

% matrix treatment
opt.mode= 'exp'; %'exp';'exact

% contrainte de temps
opt.tempsfixe_valeur = (3000)*10^-3;

% bornes
opt.Wmax   = pi/(opt.tp)/opt.mu ;
%% sert plus à rien la plus part du temps
opt.TE = 3.2e-3;
opt.TA = 0e-3 ; % ms
opt.TB =0e-3 ; %ms

opt.vec = [1];
opt.line_for_constrained_saturation = 1;
opt.spin2sature = 2; 

opt.initVec  = vecteur_initialisation(opt,samples);
opt.gradFunction = '@getGrad_ss'; 
opt.propaFunction = '@propaFunction_ss';
opt.contrainteFunction = '@getcontrainteMz';

if strcmp(opt.mode,'exact')
    opt.propaFunction = '@propafunction_ss';
    opt.initVec = convert2anglephase(opt.initVec,opt);
end

end


