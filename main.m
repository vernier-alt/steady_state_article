close all 
clear all
ref = 0;
dossier = 'segments\';
% dbstop if naninf
% DEUX ESPECES 

%% Add path

cd('C:\Users\vernier\Documents\MATLAB\steady-state\grapeSteadyState _deuxespeces_article');  
addpath(genpath('.'));



%% initialisation
[opt,samples] = set_parameters_ss(3);
opt.Instance = 0;
[spins, opt, samples] = set_grape_ss(samples,opt);
time = tic ;
Wopt = select_initialisation(spins,opt);

%% run 

opt.maxIter  = 1000; 
f0                  = 1e6;
lll = size(opt.initVec,1); % depend des paramètres a optimiser (alpha, TR, opt.Np)
it = 0;
for i=1:size(Wopt,3)
    opt.W_ini  = Wopt(:,:,i) ;
    opt.displayIter     = 'iter' ;
    [historyOpt,spins]  = grape_ss( spins, opt) ;  
    fprintf('\nCost function value : %0.4g',historyOpt.fval(end)) ;
    if (historyOpt.constraint<10^-4)  % Changer la condition
        if (historyOpt.fval(end)<f0) 
            Wopt_final       = historyOpt.W(end-lll+1:end, :) ;
            f0         =  historyOpt.fval(end)  ;
            spinsOpt   = spins ; % a voir !!!!!
            it = it + 1 ;
            c(it) = historyOpt.constraint;
        end
        
    end
    disp('iteration')
    
    i
end
opt.W_ini = Wopt_final;

%% essais 
% f0 = 1e-6;
% it = 0;
% for j = 1:numel(opt.offsetVecHz_out)
%     opt.Instance  = j;
%     [spins, opt, samples] = set_grape_ss(samples,opt);
%     opt.displayIter     = 'iter' ;
%     [historyOpt,spins]  = grape_ss(spins, opt) ;  
%     fprintf('\nCost function value : %0.4g',historyOpt.fval(end)) ;
%     if (historyOpt.constraint<10^-4)  % Changer la condition
%         if (historyOpt.fval(end)<f0) 
%             Wopt_final       = historyOpt.W(end-lll+1:end, :) ;
%             opt.W_ini        = Wopt_final;
%             f0               =  historyOpt.fval(end)  ;
%             it = it + 1 ;
%             c(it) = historyOpt.constraint;
%         end
%     end
%     disp('iteration out of')
%     j
% end
% opt.W_ini = Wopt_final;
% %%
opt.time = toc(time);


for s = 1:numel(spins)
    [spins{s}.U] = propaFunction_ss(opt.W_ini,spins{s},opt) ;
end

fprintf('\nCost function value fin initialisation : %0.4g',f0) ;

%% resultat
% opt.W_ini(1,3) = 1959e-3;
% opt.W_ini(2,3) = 19.7e-3;
% opt.W_ini(3,3) = 605.3e-3;
% 
% opt.W_ini(2,1) = 90/(opt.tp*opt.mu)*pi/180;
% opt.W_ini(3,1) = 90/(opt.tp*opt.mu)*pi/180;

result(opt);

%% affichage
h = figure;
disp('signal');
for i = 1:numel(spins)
    signal(opt.W_ini,spins{i},opt);hold on 
end
saveas(h, strcat(dossier,sprintf('%d_signal.png',ref)));
%%
Delta_f = 50; % bandwith per pixel
Delta_x = 1.172;

signal_acq = PSF(opt.W_ini,spins,opt,256,Delta_f,Delta_x,'centric',opt.W_ini,ref,dossier);
pad_excel(spins,opt,samples,opt.W_ini,signal_acq,dossier,ref)



cartographie(opt.W_ini,spins,opt,ref,dossier);

%%
load('cartos\T2_2208616.mat'); % t2 Trop haut de 7pixels
load('cartos\T1_2208616.mat');
load('cartos\DP_2208616.mat');
imshow(carto_T2,[0 max(carto_T2(:))/2]);

decalage = 1; % 7
carto_T2translate = zeros(size(carto_T2));
carto_T2translate(decalage:end,:) = carto_T2(1:end-decalage+1,:);
%imagesc(carto_T2translate*100- carto_T1);

signal_map = zeros(size(carto_T2));

for i = 1:size(carto_T2,1)
    for j = 1:size(carto_T2,2)
        spins{1}.T1 = max(carto_T1(i,j)*10^(-3), 40e-3 ); % Au cas temps de relaxation nul
        spins{1}.T2 = max(carto_T2translate(i,j)*10^(-3), 10e-3);
        spins{1}.Mt0 = 1;%carto_DP(i,j);
        [spins{1}.U] = propaFunction_ss(opt.W_ini,spins{1},opt) ;
        signal_map(i,j) = abs(get_signal(opt.W_ini,spins{1},opt,1));
    end
end
figure; 
 imshow(signal_map, [0 max(signal_map(:))*0.7])

load('cartos\T2_2016898.mat'); % t2 Trop haut de 7pixels
load('cartos\T1_2016898.mat');
load('cartos\DP_0606160.mat');