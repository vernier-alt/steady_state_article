function U = getU(optimParam,spins,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

EXPMAT  = str2func( opt.expmFunction ) ;
U = zeros(4,4,opt.Np) ;
% U(:,:,1) =  getRelaxMat(spins,optimParam(1,3),EXPMAT) * getExcMat(opt,spins,optimParam(1,1:2));
ti = opt.time_of_a_segment*optimParam(1,3)/(sum(optimParam(1:end-1,3))+opt.Nlignes*optimParam(end,3));
U(:,:,1) =  getRelaxMat(spins,ti,EXPMAT) * getExcMat(opt,spins,opt.mu*optimParam(1,1:2));



for p = 2:opt.Np
    ti = opt.time_of_a_segment*optimParam(p,3)/(sum(optimParam(1:end-1,3))+opt.Nlignes*optimParam(end,3));
    U(:,:,p) = getRelaxMat(spins,ti,EXPMAT) * getExcMat(opt,spins,opt.mu*optimParam(p,1:2))*U(:,:,p-1);
%     U(:,:,p) = getRelaxMat(spins,optimParam(p,3),EXPMAT) * getExcMat(opt,spins,optimParam(p,1:2))*U(:,:,p-1);
end



end

