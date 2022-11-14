%% Setting parameters
Ngrid = 2048;   % number of position grid points
posmax = 25;    % range of position grid
NBEC = 1e4; % number of atoms in the BEC
gBEC = 1;   % BEC intra-species interaction strength
NTGarr = 1:11; % range of particle numbers in the Tonks-Girardeau gas
gMIXarr = 0:.01:5;  % range of intra-species interaction strenghts
T = 0;  % zero temperature
[x,dx] = fftdef(posmax,Ngrid);  % defines position grid
density = 11/40; % fix density NTG/LTG for the Tonks-Girardeau gas
fname = '/data/selfpinning_zero_temperature_finite_size.mat';
fname = fullfile(fileparts(pwd),fname);

%% Changing the system size
gaps = [];
for NTG = NTGarr
    fprintf('NTG = %d:\t',NTG);
    LTG = NTG/density; % adjust size of the TG box potential

    %% Initial homogeneous BEC wave function
    wfi = sqrt(NBEC/(Ngrid*dx))*ones(Ngrid,1);

    %% Calculating the inital density for the Tonks-Girardeau gas
    %   finite difference kinetic energy operator
    TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
    wall = 1e8;
    % box potential for the TG gas
    TG_trap = zeros(Ngrid,1);
    TG_trap(abs(x)>LTG/2,:) = wall;

    % diagonalizing the single-particle Hamiltonian
    [ES,EV] = eig(TG_kin + diag(TG_trap));
    TG_energies_ini = diag(EV);
    TG_states = (1/sqrt(dx)).*ES;

    % making sure that the particle number is conserved
    Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1));
    fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; fermi_ini = fermi_ini(fermi_ini>1e-8);
    if T>0
        while abs(sum(fermi_ini) - NTG) > 1e-4
            Ef_ini = Ef_ini - T*(sum(fermi_ini) - NTG)/sum((exp((TG_energies_ini(1:length(fermi_ini)) - Ef_ini)./T).').*fermi_ini.^2);
            fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; fermi_ini = fermi_ini(fermi_ini>1e-8);
        end
    end

    % initial TG density from single-particle states
    rho_ini = sum(fermi_ini.*abs(TG_states(:,1:length(fermi_ini))).^2,2);

    clear TG_kin TG_trap ES EV TG_states

    %% Ramping up the intra-species interaction strength
    start = tic;
    for i = 1:length(gMIXarr)
        gMIX = gMIXarr(i); 
        [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini);
        % using the previous ground state as the new initial state for speed-up
        wfi = d.wf;
        rho_ini = d.rho;
        TG_energies_ini = d.TG_energies;
        params(i) = p;
        data(i) = d;
        clear p d

        gap(i) = data(i).gap(end);
    end
    toc(start)
    gaps = [gaps gap.'];
end

%% Compare with energy gap according to the analytical model 
mu0 = 0.5*gBEC*NBEC/posmax; 
for i = 1:length(gMIXarr)
    gMIX = gMIXarr(i);
    mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    a0 = 0.5*gMIX^2/gBEC;  
    a0arr(i) = a0;
    eps = 6*a0^2/(5*mu0); % we consider the limit L -> \infty resp. mubar -> mu0
    apinarr(i) = a0*(sqrt(1+2*eps)-1)/eps;
end
model = 0.5*a0arr.*apinarr;

%% Write data file for the value of the energy gap as a function of gMIX for different odd NTG in Fig. S1
header = ["gmix","NTG1","NTG3","NTG5","NTG7","NTG9","NTG11","model"];
dataout = [gMIXarr(2:end).' gaps(2:end,1:2:end) model(2:end).'];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_zero_temperature_finite_size_gap.dat'),'Delimiter','tab');

%% Write data file for the inset of Fig. S1 at fixed gMIX = 3
header = ["NTG","gap"];
dataout = [NTGarr.' gaps(301,:).'];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_zero_temperature_finite_size_gap_inset.dat'),'Delimiter','tab');

clear start i wfi rho_ini TG_energies_ini mubar a0 eps gMIX header dataout params data Ef_ini fermi_ini NTG LTG gap
save(fname);
    