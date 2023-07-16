function [Wopt] = select_initialisation(spins,opt)



%% EVALUATION D UN TRAS GRAND NOMBRE D INITIALISATIONS


nb_1 = 30;
f0_1 = 1e6*ones(nb_1,1) ; %initial cost, on en retient nb_1

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
    
    for s = 1:numel(spins)
        [spins{s}.U] = DYN(W_0,spins{s},opt) ;
    end
    pr_cost = COST(W_0,spins,opt);
    
    
    if opt.contrainte
        pr_contrainte(:) = CONT(W_0,spins,opt,1);
    else
        pr_contrainte(:) = 0;
    end

    if (max(pr_contrainte(:))<10^-2)  % if max constraint is less than 10^-2, we consider the value of the cost function
         

        if (pr_cost<max(f0_1)) % if the cost function is less that the maximum of the nb_1 previous values, we select the initialisation
            it = it + 1 ;
            num = find(f0_1==max(f0_1),1);
            Wopt_1(:,:,num)        = W_0;
            c0_1(num) = max(pr_contrainte(:));
            f0_1(num)          =  pr_cost  ;
        end

        
    elseif (max(pr_contrainte(:))<10^-1) % a second list of initializations with more tolerance is filled in case the first one would not be completed
            
       if (pr_cost<max(f0_2)) % parmis elle je retiens les 10 premières en terme de fonction de cout
            num = find(f0_2==max(f0_2),1);
            Wopt_2(:,:,num)        = W_0;
            f0_2(num)          =  pr_cost  ;
            c0_2(num) = max(pr_contrainte(:));
       end 
       
    else 
        num = find(f0_3==max(f0_3),1);
        Wopt_3(:,:,num)        = W_0;
        f0_3(num)          =  pr_cost  ;
        c0_3(num) = max(pr_contrainte(:));
    end
    i
end

[~, indices1] = sort(f0_1);
[~, indices2] = sort(f0_2);
[~, indices3] = sort(f0_3);
ll2 = size(Wopt_2,3);
if it<nb_1 % less that nb_1 initializations have been found in the first list, completion with the second one
    Wopt_1(:,:,indices1(it+1:min(nb_1,it+ll2))) = Wopt_2(:,:,indices2(1:min(ll2,nb_1-it)));
    f0_1(indices1(it+1:min(nb_1,it+ll2))) = f0_2(indices2(1:min(ll2,nb_1-it)));
    c0_1(indices1(it+1:min(nb_1,it+ll2))) = c0_2(indices2(1:min(ll2,nb_1-it)));
end
if ll2 < (nb_1-it)
    Wopt_1(:,:,it+ll2+1:nb_1) = Wopt_3(:,:,indices3(1:nb_1-ll2-it));
    f0_1(indices1(it+ll2+1:nb_1)) = f0_3(indices3(1:nb_1-ll2-it));
    c0_1(indices1(it+ll2+1:nb_1)) = c0_3(indices3(1:nb_1-ll2-it));
end

toc

f0                  = 1e6*ones(10,1) ; % initial cost, on en retient 5
opt.maxIter         = 100 ;
Wopt = zeros(size(Wopt_1,1),size(Wopt_1,2));

it = 0;
for i=1:size(Wopt_1,3)

    
    opt.W_ini  = Wopt_1(:,:,i) ;
    opt.displayIter     = 'iter' ;
    [historyOpt,spins]  = grape_ss( spins, opt) ;  
    fprintf('\n Cost function value : %0.4g',historyOpt.fval(end)) ;
    fprintf('\n Max Const Violation : %0.4g \n',historyOpt.constraint) ;
    if (historyOpt.constraint<10^-4)  % Changer la condition
        if (historyOpt.fval(end)<max(f0)) 
            Wopt(:,:,find(f0==max(f0),1))        = historyOpt.W(end-lll+1:end, :) ;
            f0(find(f0==max(f0),1))          =  historyOpt.fval(end)  ;
            it = it + 1 ;
            % c(it) = historyOpt.constraint;
        end
        
    end
%     disp('iteration')
%     i
    
end



end

