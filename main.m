close all 
clear all
ref = 0;
dossier = 'segments\';

%% Add path
% create en environnement variable for the local folder
my_path = getenv('ENV_VAR_OCMPRAGE');
if (isempty(my_path)) 
    msgbox("Create a path for the environement variable ENV_VAR_OCMPRAGE with th path to the folder (ex : C:\Documents\MATLAB\grapeSteadyState _deuxespeces_article). (Close/Reopen Matlab)");
    return;
else 
    cd(my_path);
end

addpath(genpath('.'));

%% initialization
% number of excitation pulses :
nb = 3;
[opt,samples] = set_parameters_ss(nb);
opt.Instance = 0;
[spins, opt, samples] = set_grape_ss(samples,opt);
time = tic ;
Wopt = select_initialisation(spins,opt);

%% run 
opt.maxIter  = 1000; 
f0 = 1e6;
constraint = 1e6; 
lll = size(opt.initVec,1);
for i=1:size(Wopt,3)
    opt.W_ini  = Wopt(:,:,i) ;
    opt.displayIter     = 'iter' ;
    [historyOpt,spins]  = grape_ss( spins, opt) ;  
    fprintf('\nCost function value : %0.4g',historyOpt.fval(end)) ;
    if (historyOpt.constraint<10^-4)  % Arbitrary
        if (historyOpt.fval(end)<f0) 
            Wopt_final       = historyOpt.W(end-lll+1:end, :) ;
            f0         = historyOpt.fval(end)  ;
            constraint = historyOpt.constraint;
        end 
    end
end

opt.W_ini = Wopt_final;
opt.time = toc(time);

for s = 1:numel(spins)
    [spins{s}.U] = propaFunction_ss(opt.W_ini,spins{s},opt) ;
end

fprintf('\nCost function value end : %0.4g',f0) ;
fprintf('\nCost constraint value end : %0.4g',constraint) ;

%% result
result(opt);

%% display
h = figure;
disp('signal');
for i = 1:numel(spins)
    signal(opt.W_ini,spins{i},opt);hold on 
end
saveas(h, strcat(dossier,sprintf('%d_signal.png',ref)));
%%
Delta_f = 50; % bandwith per pixel
Delta_x = 1.172;

signal_acq = PSF(opt.W_ini,spins,opt,opt.Nlignes*opt.Ncylces,Delta_f,Delta_x,'centric',opt.W_ini,ref,dossier);
pad_excel(spins,opt,samples,opt.W_ini,signal_acq,dossier,ref)
cartographie(opt.W_ini,spins,opt,ref,dossier);

