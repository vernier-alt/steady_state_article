function [] = cartographie(W,spins,opt,ref,dossier)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
DYN             = str2func( opt.propaFunction ) ;

T1tab = [];
T2tab = [];
for i = 1:numel(spins)
    T1tab = horzcat(T1tab, spins{i}.T1);
    T2tab = horzcat(T2tab, spins{i}.T2);
end 


couleur  = ["#7E2F8E"; "#EDB120" ;"#77AC30";"#D95319"];

% graisse,

valeurs_T1 = linspace(0.05,1.2*max(T1tab),50);
valeurs_T2 = linspace(0.005,1.2*max(T2tab),50);
% valeurs_T1 = linspace(0.05,2000e-3,50);
% valeurs_T2 = linspace(0.005,200e-3,50);
% 
% valeurs_T1 = linspace(0.05,2500e-3,50);
% valeurs_T2 = linspace(0.005,150e-3,50);
% 
% valeurs_T1 = linspace(0.05,4500e-3,50);
% valeurs_T2 = linspace(0.005,2500e-3,50);

out = zeros(numel(valeurs_T2),numel(valeurs_T1));

spins_memoire = spins;
opt_Niso_memory = opt.N_iso;

    for k=1:numel(valeurs_T2)
        for i=1:numel(valeurs_T1)

            if opt.time_of_a_segment < 0
                disp('error');
            end 
            samples.T2 = [valeurs_T2(k)];
            samples.T1 = [valeurs_T1(i)];
            samples.PD = [1];
            samples.maxi = [1];
            samples.saturation = [0];

            
            [spins, opt, samples] = set_grape_ss(samples,opt);


            for s = 1:numel(spins)
                [spins{s}.U] = DYN(W,spins{s},opt) ;
            end

%             out(k,i) = (abs(get_signal(spins{1},opt)) -abs(get_signal(spins{2},opt)))/(abs(get_signal(spins{1},opt)) +abs(get_signal(spins{2},opt)));
%               for ligne = 1:opt.Nlignes
%                   
                ligne = 1;
                out(k,i) = out(k,i) + get_signal(W,spins{1},opt,ligne);
%               end 
              out(k,i) = abs(out(k,i));

        end 
    end 

%%
% l = {};
% clear a 
% clear legh
% for j=1:numel(segment_duration)
%     a(j)=plot(valeurs_T1,out(j,:)); hold on 
%     legh{j}=sprintf('segment duration %s ms',num2str(segment_duration(j)*1000));
% end 
% legend(a,legh);
% xlabel('T1 ms, reference 1945 ms');
%%

% mini = min(out(:));
 mini = 0;
% maxi = max(out(:));
 maxi = 0.2;
step = 10;
ncouleurs = round((maxi - mini)/0.001);%0.001
% ncouleurs = 100;

%% 

k = figure;

    colormap(jet(ncouleurs))

    [T1 T2] = meshgrid(valeurs_T1*10^3, valeurs_T2*10^3);
    contourf(T1,T2,out(:,:),30); hold on
    xlabel('T1 ms');
    ylabel('T2 ms');
    hold on

    for num = 1:opt_Niso_memory

            scatter(spins_memoire{num}.T1*10^3,spins_memoire{num}.T2*10^3,100,'filled','MarkerFaceColor',couleur(num));
            hold on
%          end
    end 
    caxis([mini maxi])
    colorbar
    
%saveas(k, strcat('manipes_2704\',sprintf('%d_carto.png',ref)));    
saveas(k, strcat(dossier,sprintf('%d_carto.png',ref)));  
  
    
end









