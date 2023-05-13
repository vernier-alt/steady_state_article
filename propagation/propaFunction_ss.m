function [U] = propaFunction_ss(optimParam,spins,opt)

% disp('Bonjourrrrrrrr')


EXPMAT  = str2func( opt.expmFunction ) ;

B0      = 0;


U = zeros(4,4,opt.Np+2);
U(:,:,1) = eye(4);

if strcmp(opt.mode,'exp')
    optimParam(1:end,1:2) = opt.mu*optimParam(1:end,1:2);
else
    optimParam(1:end,1) = opt.mu*optimParam(1:end,1);
end 

ti = getDelais_changementvar(opt,optimParam);


for p = 2:opt.Np+1
   U(:,:,p) = getRelaxMat(spins,ti(p-1),EXPMAT,opt.mode) * getExcMat(opt,spins,optimParam(p-1,1:2))*U(:,:,p-1);
end

U = U(:,:,2:opt.Np+1);




end 

