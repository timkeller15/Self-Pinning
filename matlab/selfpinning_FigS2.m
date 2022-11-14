%% Setting parameters
Ngrid = 2048;   % number of position grid points
posmax = 25;    % range of position grid
density = 11/40; % fix density NTG/LTG for the Tonks-Girardeau gas
NBEC = 1e4; % number of atoms in the BEC
gBEC = 1;   % BEC intra-species interaction strength
gMIX = 3;  % inter-species interaction strenght
temperatures = 0:.01:5; % temperature range
NTGarr = 1:2:11; % range of system sizes
[x,dx] = fftdef(posmax,Ngrid);  % defines position grid
wall = 1e8; % height of box potential 

pinnedness_arr = [];
mompeaks_arr = [];

%% Changing the system size
% The outer loop over the system sizes should be parallelized as
% much as possible for reasonable execution times. For example by commenting it out
% and setting NTG = NTGarr(${SLURM_ARRAY_TASK_ID}) when running the code
% as an array job on a cluster using slurm. 

for NTG = NTGarr
    fprintf('NTG = %d:\n',NTG);
    LTG = NTG/density; % adjust size of the TG box potential
    
    fname = sprintf('/data/selfpinning_finite_temperature_finite_size_NTG%d_gMIX3_momentum.mat',NTG);
    fname = fullfile(fileparts(pwd),fname);
    
    % resetting inital parameters 
    TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
    TG_trap = zeros(Ngrid,1);
    TG_trap(abs(x)>LTG/2,:) = wall;

    [ES,EV] = eig(TG_kin + diag(TG_trap));
    TG_energies_ini = diag(EV);
    TG_states = (1/sqrt(dx)).*ES;

    T = 0;
    Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1));
    fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; fermi_ini = fermi_ini(fermi_ini>1e-8);
    rho_ini = sum(fermi_ini.*abs(TG_states(:,1:length(fermi_ini))).^2,2);
    clear Ef_ini fermi_ini

    wfi = sqrt(NBEC/(Ngrid*dx))*ones(Ngrid,1);

    clear TG_kin TG_trap ES EV TG_states

    momdis = zeros(Ngrid,length(temperatures));
    k = (pi/dx)*linspace(0,1,Ngrid/2 + 1); 
    k = [-fliplr(k) k(2:end-1)]; 
    dk = pi/posmax;

    %% Calculate the momentum distribution from the reduced single-particle density matrix
    % Only the lower triangle of the RSPDM is actually computed and then mirrored, since the matrix is symmetric  

    for i = 1:length(temperatures) 
        T = temperatures(i); 
        fprintf('T = %3.2f:\t',T);
        [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini);
        wfi = d.wf;
        rho_ini = d.rho;
        TG_energies_ini = d.TG_energies;
        params(i) = p;
        data(i) = d;
        clear p d
        
        pinnedness_num(i) = sum(data(i).fermi(1:NTG))/NTG;

        tic; rspdm = selfpinning_rspdm(Ngrid,posmax,LTG,gMIX,temperatures(i),params(i).Ef,data(i).wf); toc
        [ES,EV] = eig(rspdm);
        occupationnr = flipud(diag(EV))*dx;
        naturalorbitals=(1/sqrt(dx))*fliplr(ES);
        occupationnr = occupationnr(occupationnr>0);
        naturalorbitals = naturalorbitals(:,occupationnr>0);
        momdis(:,i) = fftshift(abs(fft(naturalorbitals)*(dx/sqrt(2*pi))).^2*occupationnr); % normalisation factor / Parseval's theorem 
        clear ES EV rspdm occupationnr naturalorbitals
        save(fname, '-regexp', '^(?!(pinnedness_arr|mompeaks_arr)$).');
    end
    clear i wfi rho_ini TG_energies_ini T  
    save(fname, '-regexp', '^(?!(pinnedness_arr|mompeaks_arr)$).');
    
    pinnedness_arr = [pinnedness_arr pinnedness_num.'];
    mompeaks_arr = [mompeaks_arr momdis(Ngrid/2 + 1,:).'];
    
    clear params data momdis pinnedness_num 
end

a0 = 0.5*gMIX^2/gBEC; mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax; eps = 6*a0^2/(5*mubar); ap = a0*(sqrt(1+2*eps)-1)/eps;
Tf = 0.5*a0*ap; % scale temperature in units of the energy gap size at T = 0 
%% Write data file for the ground state occupancy in Fig. S2 (a)
header = ["T","NTG1","NTG3","NTG5","NTG7","NTG9","NTG11"];
dataout = [temperatures.'/Tf pinnedness_arr];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_finite_size_pinnedness.dat'),'Delimiter','tab');

%% Write data file for the zero-momentum peaks in Fig. S2 (b) 
header = ["T","NTG1","NTG3","NTG5","NTG7","NTG9","NTG11"];
dataout = [temperatures.'/Tf mompeaks_arr];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_finite_size_momentum_peaks.dat'),'Delimiter','tab');

clear a0 mubar eps ap header dataout
fname = sprintf('/data/selfpinning_finite_temperature_finite_size_gMIX3_momentum.mat',NTG);
fname = fullfile(fileparts(pwd),fname);
save(fname);