function [] = result(opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ordre={};
for i=1:opt.Np
    ordre=union(ordre,sprintf('%d',i));
end

if opt.changement_of_variable
    Delay_ms=opt.W_ini(1:opt.Np,3)*10^3*(opt.time_of_a_segment-opt.TR*opt.Nlignes)/(sum(opt.W_ini(1:opt.Np,3)));
    TR = opt.TR*10^3;
else
    Delay_ms=opt.W_ini(1:opt.Np,3)*10^3;
    TR = opt.TR*10^3; 
end

if strcmp(opt.mode,'exp')
    disp('\n début');
    RFx=floor(opt.W_ini(1:opt.Np,1)*opt.tp*180/pi*10*opt.mu)/10;
    RFy=floor(opt.W_ini(1:opt.Np,2)*opt.tp*180/pi*10*opt.mu)/10;
    init = table(ordre', RFx,RFy,Delay_ms(1:opt.Np))
    fprintf('L angle optimal est %0.2f degres et le TR optimal est de %0.1f ms \n',opt.alpha*180/pi,TR);

elseif strcmp(opt.mode,'exact')
    disp('\n début');
    RF=floor(opt.W_ini(1:opt.Np,1)*opt.tp*180/pi*10*opt.mu)/10;
    Phi=floor(opt.W_ini(1:opt.Np,2)*180/pi*10)/10;
    init = table(ordre', RF,Phi,Delay_ms(1:opt.Np))
    fprintf('L angle optimal est %0.2f degres et le TR optimal est de %0.1f ms \n',opt.alpha*180/pi,TR);


end

end

