function [Wopt_2] = select_initialisation(spins,opt)



%% EVALUATION D UN TRAS GRAND NOMBRE D INITIALISATIONS


nb_1 = 30;
f0_1 = 1e6*ones(nb_1,1) ; %initial cost, on en retient 10

f0_2 = 1e6*ones(nb_1,1) ; 
f0_3 = 1e6*ones(nb_1,1) ; 


lll = size(opt.initVec,1); % depend des paramètres a optimiser (alpha, TR, opt.Np)
Wopt_1 = zeros(size(opt.initVec,1),size(opt.initVec,2));
Wopt_2 = Wopt_1;
Wopt_3 = Wopt_1;

c0_1 = f0_1;
c0_2 = f0_2;
c0_3 = f0_3;

COST    = str2func( opt.costFunction );
DYN = str2func( opt.propaFunction );
CONT = str2func( opt.contrainteFunction );

it = 0;

tic

for i=1:size(opt.initVec,3)

    W_0  = opt.initVec(:,:,i) ;
    
    % pour aller plus vite, pré-selection à la resonnance
    for s = 1:numel(spins)
        [spins{s}.U] = DYN(W_0,spins{s},opt) ;
    end
    pr_cost = COST(W_0,spins,opt);
    
    
    if opt.contrainte
        pr_contrainte(i,:) = CONT(W_0,spins,opt,1);
    else
        pr_contrainte(i,:) = 0;
    end

    if (max(pr_contrainte(i,:))<10^-2)  % Si la containte est acceptable à 10^-2 près, on garde
         
        test = true; 
%         if opt.optimisation_zcoherence % si la deuxième contrainte est activée et que l'on la repecte on garde 
%             for s = 1:numel(spins)
%                 [spins{s}.U] = DYN(W_0,spins{s},opt) ;
%             end
%             test = (max(z_coherence(W_0,spins,opt))<0.1); 
%         end
        if test
            if (pr_cost<max(f0_1)) % parmis elle je retiens les 10 premières en terme de fonction de cout
                it = it + 1 ;
                num = find(f0_1==max(f0_1),1);
                Wopt_1(:,:,num)        = W_0;
                c0_1(num) = max(pr_contrainte(i,:));
                f0_1(num)          =  pr_cost  ;
                
            end
            
        end
        
    elseif (max(pr_contrainte(i,:))<10^-1)
            
       if (pr_cost<max(f0_2)) % parmis elle je retiens les 10 premières en terme de fonction de cout
            num = find(f0_2==max(f0_2),1);
            Wopt_2(:,:,num)        = W_0;
            f0_2(num)          =  pr_cost  ;
            c0_2(num) = max(pr_contrainte(i,:));
       end 
       
    else 
        num = find(f0_3==max(f0_3),1);
        Wopt_3(:,:,num)        = W_0;
        f0_3(num)          =  pr_cost  ;
        c0_3(num) = max(pr_contrainte(i,:));
    end
    i
end

[~, indices2] = sort(f0_2);
[~, indices3] = sort(f0_3);
ll2 = size(Wopt_2,3);
if it<nb_1
    Wopt_1(:,:,it+1:min(nb_1,it+ll2)) = Wopt_2(:,:,indices2(1:min(ll2,nb_1-it)));
    f0_1(it+1:min(nb_1,it+ll2)) = f0_2(indices2(1:min(ll2,nb_1-it)));
    c0_1(it+1:min(nb_1,it+ll2)) = c0_2(indices2(1:min(ll2,nb_1-it)));
end
if ll2 < (nb_1-it)
    Wopt_1(:,:,it+ll2+1:nb_1) = Wopt_3(:,:,indices3(1:nb_1-ll2-it));
    f0_1(it+ll2+1:nb_1) = f0_3(indices3(1:nb_1-ll2-it));
    c0_1(it+ll2+1:nb_1) = c0_3(indices3(1:nb_1-ll2-it));
end

toc

f0                  = 1e6*ones(10,1) ; % initial cost, on en retient 5
opt.maxIter         = 100 ;
Wopt_2 = zeros(size(opt.initVec,1),size(opt.initVec,2));

it = 0;
for i=1:size(Wopt_1,3)

    
    opt.W_ini  = Wopt_1(:,:,i) ;
    opt.displayIter     = 'iter' ;
    [historyOpt,spins]  = grape_ss( spins, opt) ;  
    fprintf('\n Cost function value : %0.4g',historyOpt.fval(end)) ;
    fprintf('\n Max Const Violation : %0.4g \n',historyOpt.constraint) ;
    if (historyOpt.constraint<10^-4)  % Changer la condition
        if (historyOpt.fval(end)<max(f0)) 
            Wopt_2(:,:,min(find(f0==max(f0))))        = historyOpt.W(end-lll+1:end, :) ;
            f0(min(find(f0==max(f0))))          =  historyOpt.fval(end)  ;
%             spinsOpt(min(find(f0==max(f0))))   = spins ; % a voir !!!!!
            it = it + 1 ;
            c(it) = historyOpt.constraint;
        end
        
    end
%     disp('iteration')
%     i
    
end



end

