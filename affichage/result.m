function [] = result(opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ordre={};
for i=1:opt.Np
    ordre=union(ordre,sprintf('%d',i));
end

if opt.optimisation_TR
    if opt.changement_variable_delais
        Delay_ms=opt.W_ini(1:end,3)*10^3*opt.tempsfixe_valeur/(sum(opt.W_ini(1:end-1,3))+opt.Nlignes*opt.W_ini(end,3));
        TR = Delay_ms(end);
    else
        Delay_ms=opt.W_ini(1:end,3)*10^3;
        TR = Delay_ms(end);
    end
else
    if opt.changement_variable_delais
        Delay_ms=opt.W_ini(1:opt.Np,3)*10^3*(opt.tempsfixe_valeur-opt.TR*opt.Nlignes)/(sum(opt.W_ini(1:opt.Np,3)));
        TR = opt.TR*10^3;
    else
        Delay_ms=opt.W_ini(1:opt.Np,3)*10^3;
        TR = opt.TR*10^3; 
    end
end
if strcmp(opt.mode,'exp')
    disp('\n début');
    RFx=floor(opt.W_ini(1:opt.Np,1)*opt.tp*180/pi*10*opt.mu)/10;
    RFy=floor(opt.W_ini(1:opt.Np,2)*opt.tp*180/pi*10*opt.mu)/10;
    init = table(ordre', RFx,RFy,Delay_ms(1:opt.Np))
    fprintf('L angle optimal est %0.2f degres et le TR optimal est de %0.1f ms \n',getalpha_changementvar(opt,opt.W_ini)*180/pi,TR);

elseif strcmp(opt.mode,'exact')
    disp('\n début');
    RF=floor(opt.W_ini(1:opt.Np,1)*opt.tp*180/pi*10*opt.mu)/10;
    Phi=floor(opt.W_ini(1:opt.Np,2)*180/pi*10)/10;
    init = table(ordre', RF,Phi,Delay_ms(1:opt.Np))
    fprintf('L angle optimal est %0.2f degres et le TR optimal est de %0.1f ms \n',getalpha_changementvar(opt,opt.W_ini)*180/pi,TR);


end

end

