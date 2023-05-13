function [signal_acq] = PSF(optimParam,spins,opt,n_lignes_matrice,delta_f,delta_x,mode,W,ref,dossier)

alpha = getalpha_changementvar(opt,W);
TR = getTR_changementvar(opt,W);


TR = TR*10^3; %ms 
TA = opt.TA *10^3; % ms
TB = opt.TB *10^3; %ms

% if strcmp(mode,'centric')
%     rapport = n_lignes_matrice/opt.Nlignes;
%     centre_kspace = floor(
%     order = [
% end 
rapport = n_lignes_matrice/opt.Nlignes;
signal_acq = ones(numel(spins),n_lignes_matrice,1);
centre_kspace = floor(n_lignes_matrice/2);
x_espace_direct = [-floor(n_lignes_matrice/2):1:floor(n_lignes_matrice/2)-1]*delta_x;
x_espace_indirect = [-floor(n_lignes_matrice/2):1:floor(n_lignes_matrice/2)-1];

n1 = 2^nextpow2(length(x_espace_direct));

k = figure;

for num=1:numel(spins)
    
    [spins{num}.U]   = propaFunction_ss(optimParam,spins{num},opt) ; 

    T1 = spins{num}.T1*10^3; %ms
    T2 = spins{num}.T2*10^3 ; % ms 

    n = opt.Nlignes ;
    n_cycle = opt.Ncylces;



    H_U = [0 0 0 1]*spins{num}.U(:,:,opt.Np)*[1 0 0 0]';
    F_U =  [0 0 0 1]*spins{num}.U(:,:,opt.Np)*[0 0 0 1]';

    % EXP
    EA = exp(-TA/T1); %ms
    EB = exp(-TB/T1); %ms
    E1 = exp(-TR/T1); %ms
    K = cos(alpha)*E1; %ms 
    KK =  K^n*EA*EB*F_U;
    
%     nbinter = 1;
% time_TR = [0:floor(TR/nbinter)-1];

    M0 = spins{num}.Mt0(end);
    % STEADY - STATE 
    ss = (M0*((1-EB)+( (1-E1)*(1-K^n)/(1-K)+K^n*(1-EA)))*EB + K^n*EA*EB*H_U)/(1-K^n*EA*EB*F_U);
    s = M0*(1-E1)/(1-K);
    ss_EB = (ss -(1-EB))/EB;


    Mz = s +(ss_EB-s)*K^(-opt.Nlignes);
    
    
    
    for j=1:opt.Nlignes
        kk_haut = centre_kspace - (j-1)*rapport/2;
        kk_bas = centre_kspace + (j-1)*rapport/2 +1 ;
%         signal_acq(num,kk_haut - (rapport/2 -1):kk_haut ) =  (M0 + ( Mz*cos(alpha) - M0 )*exp(-TR/T1))*sin(alpha)*exp(-opt.TE/spins{num}.T1);
%         signal_acq(num,kk_bas:kk_bas + (rapport/2 -1)) =  (M0 + ( Mz*cos(alpha) - M0 )*exp(-TR/T1))*sin(alpha)*exp(-opt.TE/spins{num}.T1);
        signal_acq(num,kk_haut - (rapport/2 -1):kk_haut ) =  Mz*sin(alpha)*exp(-opt.TE/spins{num}.T1);
        signal_acq(num,kk_bas:kk_bas + (rapport/2 -1)) =  Mz*sin(alpha)*exp(-opt.TE/spins{num}.T1);
        Mz = M0 + ( Mz*cos(alpha) - M0 )*exp(-TR./T1);
    end
    
    subplot(1,2,1);
    plot(x_espace_indirect,signal_acq(num,:));hold on
    xlabel(' Lines in k space');
    subplot(1,2,2);
    plot(x_espace_direct,ifftshift(real(ifft(signal_acq(num,:)')))); hold on 
%     plot(x_espace_direct,ifftshift((ifft(signal_acq(num,:)')))); hold on 
    r1 = sum(ifftshift(abs(ifft(signal_acq(num,:)'))))/(numel(signal_acq(num,:)));
    r2 = sum(ifftshift((ifft(signal_acq(num,:)'))))/(numel(signal_acq(num,:)));
    xlim([-10 10]);
    ylim([-0.05 0.15]);
    xlabel(' Y (mm)');
    
    integrale = sum(signal_acq(num,:))
end

saveas(k, strcat(dossier,sprintf('%d_PSF.png',ref)));

% figure;
% plot(ifftshift(abs(ifft(signal_acq(1,:)')))-ifftshift(abs(ifft(signal_acq(2,:)'))));

end

