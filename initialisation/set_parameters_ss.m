function [opt,samples] = set_parameters_ss(k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% tubes 
% tissue1 = [1500 380 1300 3500]*10^-3;
% tissue2 = [80 140 35 350]*10^-3;
% tissue0 = [1 1 1 1];

tubes1 = [120 106    65    86    139    74    483    216    633    359    762   770]*10^-3;
tubes2 = [9   11     8.5   9       8    10     32     16     55     27     72     46]*10^-3;
tubes0 = [1   0.9    0.8   0.7   1      1      1      1      1      1        1      1 ];
muscleT2 = 34e-3;
muscleT1 = 1249e-3;
fibroseT2 = 64e-3;
fibroseT1 = 1912e-3 ;
tissussainT1 = 1994e-3;
tissussainT2 = 94e-3;

tumeurT1 = 2100e-3;
tumeurT2 = 400e-3;

% T1,T2
samples.T1 = [tubes1(11) tubes1(7) tubes1(6)];
samples.T2 = [tubes2(11) tubes2(7) tubes2(6)];
samples.PD = [tubes0(11) tubes0(7) tubes0(6)];


samples.T1 = [tissussainT1 tumeurT1  ];
samples.T2 = [tissussainT2 tumeurT2 ];
samples.PD = [1 1];


% % 
samples.maxi = [1 0  0];
samples.saturation = [0 1 1];

% samples.T1 = [tubes1(9) tubes1(11) ];
% samples.T2 = [tubes2(9) tubes2(11) ];
% 
% samples.PD = [tubes0(9) tubes0(11) ];

% GM , WM, CSF
brainT1 = [1450 900 4200]*10^-3;
brainT2 = [96 70 2000]*10^-3;
brainM0 = [0.75 0.65 1];

samples.T1 = [brainT1(1) brainT1(2) brainT1(3)];
samples.T2 = [brainT2(1) brainT2(2) brainT2(3)];

samples.PD = [brainM0(1) brainM0(2) brainM0(3)];

 samples.maxi = [1 0 0];
 samples.saturation = [0 1 1];

% T1 = 74.6;
% T2 = 10;
% 
% samples.T1 = [T1*15/10       T1]*10^-3;
% samples.T2 = [T2*15/10  T2]*10^-3;
% 
% samples.PD = [1 1];
% 
% samples.maxi = [0 1];
% samples.saturation = [1 0];
% samples.PD = [0.75 0.65];
% samples.T1= [460 450]*10^-3;%cc, cx
% samples.T2 = [8 18]*10^-3;
% samples.maxi = [1 0];
% samples.saturation = [0 1];


opt.saturation_tolerance = 0.0001; % tolerance 0.001 c'est bien


% opt.nbFreq         = 31 ;
% opt.BandWidthMax      = 100 ; % HZ
% pd = makedist('Normal','mu',0,'sigma',opt.BandWidthMax/4); % crée une distribution gaussienne 
% u = cdf(pd,opt.BandWidthMax/2); % calcule la cdf pour atteindre BWMAX 
% x = linspace(0.5,u,opt.nbFreq); % crée un echantillonnage
% opt.offsetVecHz = icdf(pd,x); %inverse
opt.offsetVecHz_out = [10 12 25 50 75 125 130 199 253 325 409 450];%[199 253 409];%[199 253 409]; % 59 1507 809];


opt.Np = k+1;
opt.tp     = 1e-3;%1e-3; % pulse duration (in ms)
opt.expmFunction        = '@expmdemo1';
opt.costFunction = '@getWG_epsilon';%'@getWG_epsilon_B0''@getWG_epsilon'; % '@getCO', '@getWGABS';'@get_signal_carre'; %'@get_ss_carre'; %'@get_signal_carre'

opt.sat = 1; % pour getCO
opt.epsilon_abs = 10^-5; % pour wg_epsilon

%% pour Grape
% opt.TolX                = 1e-8 ; % Step Tolerance
% opt.TolFun              = 1e-8 ; % Function Tolerance

opt.optimisation_zcoherence = false;

opt.TR_alpha = false ;
opt.epsilon_alpha_TR     = 1e-7 ;

opt.TolX                = 1e-15 ; % Step Tolerance
opt.TolFun              = 1e-5 ; % Function Tolerance
opt.epsilon     = 1e-10 ;

opt.displayIter = 'iter' ;
opt.optimFunction       = 'fmincon';% 'fmincon' ; 

% Les cahngements de variables 
opt.changement_variable_delais = true ; % si changement de variable, on conserve une contrainte d'inégalité si TR min
opt.mu = 1000; % b = mu * optparam, si mu trop grand, optparam trop petit, gradient élevé ereur de derivative check

% Les paramètres d'acquisition
opt.optimisation_TR = false;
opt.TR = 6e-3; % Sinon voici le TR

opt.optimisation_alpha = false;
opt.alpha = 12*pi/180;% Sinon voici le alpha 

opt.Kmin = cos(12*pi/180)*exp(-6e-3/1.5);

% matrix treatment
opt.mode= 'exp'; %'exp';'exact

% contrainte de temps
opt.tempsfixe_valeur = (5000)*10^-3;
% Lignes
opt.Nlignes = 64;
opt.Ncylces = 30;% graphique

% bornes
opt.Wmax   = pi/(opt.tp)/opt.mu ;

opt.TRmax = 12e-3;
opt.TRmin = 6e-3;

opt.alphamax = 15*pi/180;
opt.alphamin = 10*pi/180 ;

%% sert plus à rien la plus part du temps
opt.TE = 3.2e-3;
opt.TA = 0e-3 ; % ms
opt.TB =0e-3 ; %ms


opt.vec = [1];
opt.spin2sature = 2;
opt.target_z = 0;



opt.initVec  = vecteur_initialisation(opt,samples);

% if strcmp(opt.mode,'exact')
%     opt.propaFunction = '@propaFunction_ss_exact'; % '@propaFunction_ss'; % 'propafunction_ss_exact'
%     opt.gradFunction = '@getGrad_ss_exact'; % '@GetGrad_ss'
%     opt.gradFunction = '@getGrad_ss_exact'; % '@GetGrad_ss'
%     opt.contrainteFunction = '@getcontrainteMz';
% else
    opt.propaFunction = '@propaFunction_ss'; % '@propaFunction_ss'; % 'propafunction_ss_exact'
    opt.gradFunction = '@getGrad_ss'; % '@GetGrad_ss'
    opt.contrainteFunction = '@getcontrainteMz';
% end 

if strcmp(opt.mode,'exact')
    opt.initVec = convert2anglephase(opt.initVec,opt);
end

end


