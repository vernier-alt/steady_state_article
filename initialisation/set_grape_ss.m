function [spins, opt, samples]= set_grape_ss(samples,opt)

%% Set_grape_lagrange(samples);
opt.N_iso       = max([numel(samples.T1),numel(samples.T2),numel(samples.PD)]);
s = 1; 

opt.contrainte = (max(samples.saturation)>0) ;

for i = 1:opt.N_iso 
    
        iso.Mt0         = [1 ; 0 ; 0 ; samples.PD(i)] ; % equilibrium magnetization state [1 ; Mx ; My ; Mz] depending on proton density

        iso.T       	= [0 0 0 0 ; 0 -1/samples.T2(i) 0 0 ; 0 0 -1/samples.T2(i) 0 ; samples.PD(i)/samples.T1(i) 0 0 -1/samples.T1(i)]; %relaxation terms

        iso.B1_inh = 1; % for further implementations, B1 inhomogeneities
        iso.B0_inh = 0; % for further implementations, B0 inhomogeities

        iso.Wx      = iso.B1_inh*[0 0 0 0 ; 0 0 0 0 ; 0 0 0 +1 ; 0 0 -1 0] ; %W1x matrix
        iso.Wy      = iso.B1_inh*[0 0 0 0 ; 0 0 0 -1 ; 0 0 0 0 ; 0 +1 0 0] ; %W1y matrix

        iso.T1       	= samples.T1(i) ;
        iso.T2       	= samples.T2(i) ;

        iso.W0          = [0 0 0 0 ; 0 0 +1 0 ; 0 -1 0 0 ; 0 0 0 0] ; %The W0 matrix for propagation
        
        if samples.maxi(i)==0  
            iso.max= -1; % saturer
        else
            iso.max = 1; % max
        end

        if (samples.saturation(i)== 1)
            iso.saturation = true; % saturation states
        else
            iso.saturation = false;
        end
        
        spins{s} = iso; 
        s = s + 1 ; %increment spin index
    end
end
