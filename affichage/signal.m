function [] = signal(optimParam,spins,opt)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
% données
TR = getTR_changementvar(opt,optimParam)*10^3;
T1 = spins.T1*10^3; %ms
% T2 = spins.T2*10^3 ; % ms 

TA = opt.TA *10^3; % ms
TB = opt.TB *10^3; %ms

n = opt.Nlignes ;
n_cycle = opt.Ncylces;

alpha = getalpha_changementvar(opt,optimParam);

H_U = [0 0 0 1]*spins.U(:,:,opt.Np)*[1 0 0 0]';
F_U =  [0 0 0 1]*spins.U(:,:,opt.Np)*[0 0 0 1]';

% EXP
EA = exp(-TA/T1); %ms
EB = exp(-TB/T1); %ms
E1 = exp(-TR/T1); %ms
K = cos(alpha)*E1; %ms 
KK =  K^n*EA*EB*F_U;



% STEADY - STATE 
ss = (spins.Mt0(end)*((1-EB)+( (1-E1)*(1-K^n)/(1-K)+K^n*(1-EA)))*EB+K^n*EA*EB*H_U)/(1-K^n*EA*EB*F_U);


s = spins.Mt0(end)*(1-E1)/(1-K);
ss_EB = (ss -(1-EB))/EB;

if (spins.B0_inh == 0)
    signal = (s + (ss_EB - s ) * K^(-floor(opt.Nlignes)))*exp(-opt.TE/spins.T2)*sin(alpha)
else
    signal = (s + (ss_EB - s ) * K^(-floor(opt.Nlignes)))*exp(-opt.TE/spins.T2)*sin(alpha);
end

nbinter = 1; % nombre de points par ms 
% temps_cycle = TA*n*TR+TB;
time_TA = [0:floor(TA/nbinter)-1];
time_TB = [0:floor(TB/nbinter)-1];
time_TR = [0:floor(TR/nbinter)-1];

S = [1];
M0 = spins.Mt0(end);
Mz = M0;
signal = [M0];



for i=1:n_cycle
    signal = horzcat(signal,M0 + ( F_U*Mz+ H_U - M0 )*exp(-time_TA./T1));
    Mz = M0 + ( F_U*Mz+ H_U - M0 )*exp(-TA./T1);
    for j=1:n
        signal = horzcat(signal,M0 + ( Mz*cos(alpha) - M0 )*exp(-time_TR./T1));
        Mz = M0 + ( Mz*cos(alpha) - M0 )*exp(-TR./T1);
    end
    signal = horzcat(signal,M0 + ( Mz - M0 )*exp(-time_TB./T1));
    Mz = M0 + ( Mz - M0 )*exp(-TB./T1);
end

% conv = - EA*EB*(cos(alpha)*E1)^n;

temps = [0 : size(signal,2)-1]/nbinter;
enveloppe1 = exp((temps)/(TA+n*TR+TB)*log(abs(KK)))*(M0-ss)+ss;
enveloppe2 = -exp((temps)/(TA+n*TR+TB)*log(abs(KK)))*(M0-ss)+ss;

plot(0*signal); hold on 
plot(signal); hold on
plot(ss*ones(size(signal)));hold on 
plot(enveloppe1);
plot(enveloppe2);

legend('0','signal','steady-state','enveloppemax','eveloppemin');

end

