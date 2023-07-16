function [] = pad_excel(spins,opt,samples,W,signal_acq,dossier,ref)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

alpha = opt.alpha
TR = opt.TR;

A ={};

A = horzcat(A,datestr(datetime));
A = horzcat(A,{' '});
A = horzcat(A,ref);
A = horzcat(A,{' '});
A = horzcat(A,opt.mode);
A = horzcat(A,{' '});
A = horzcat(A,opt.time);
A = horzcat(A,{' '});

for j=1:4
    if j<=numel(spins)
        A = horzcat(A,spins{j}.T1*10^3);
        A = horzcat(A,spins{j}.T2*10^3);
    else
        A = horzcat(A,{' '});
    end
end 

A = horzcat(A,{' '});
A = horzcat(A,opt.Np-1);
A = horzcat(A,opt.time_of_a_segment*10^3);
A = horzcat(A,opt.Nlignes);
A = horzcat(A,TR*10^3);
A = horzcat(A,alpha*180/pi);
A = horzcat(A,opt.contrainte);
A = horzcat(A,{' '});


del = getDelais_changementvar(opt,W);

A = horzcat(A,del(1)*10^3);

for i = 2:5
    if i<=opt.Np
        if strcmp(opt.mode,'exact')
            A = horzcat(A,floor(W(i,1)*opt.tp*180/pi*10*opt.mu)/10);
            A = horzcat(A,floor(W(i,2)*180/pi*10)/10);
        else
           A = horzcat(A,floor(sqrt(W(i,1)^2+W(i,2)^2)*opt.tp*180/pi*10*opt.mu)/10);
           A = horzcat(A,floor(angle(W(i,1)+1i*W(i,2))*10)/10);
        end
        A = horzcat(A,del(i)*10^3);
    else
        A = horzcat(A,{' '});
    end
end 

A = horzcat(A,{' '});

for k = 1:numel(spins)
    A = horzcat(A,get_signal(W,spins{k},opt,1));
end 

 A = horzcat(A,{' '});
 
for k = 1:numel(spins)
    A = horzcat(A,sum(signal_acq(k,:)));
end 

fichier = strcat(dossier,'donnees.xlsx');
xlsappend(fichier,A);

% xlsappend('manipes_2704\essais_2704.xlsx',A);

end


