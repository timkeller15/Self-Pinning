%% Setting parameters
Ngrid = 2048;   % number of position grid points
posmax = 25;    % range of position grid
NTG = 7;    % number of atoms in the Tonks-Girardeau gas
LTG = 40;   % length of the box potential for the TG gas
NBEC = 1e4; % number of atoms in the BEC
gBEC = 1;   % BEC intra-species interaction strength
gMIXarr = 0:.01:5;  % range of inter-species interaction strenghts
T = 0;  % zero temperature
[x,dx] = fftdef(posmax,Ngrid);  % defines position grid
fname = '/data/selfpinning_zero_temperature_NTG7_gMIXramp.mat';
fname = fullfile(fileparts(pwd),fname);

%% Initial homogeneous BEC wave function
wfi = sqrt(NBEC/(Ngrid*dx))*ones(Ngrid,1);

%% Calculating the inital density for the Tonks-Girardeau gas
% finite difference kinetic energy operator
TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));

% box potential for the TG gas
TG_trap = zeros(Ngrid,1);
wall = 1e8;
TG_trap(abs(x)>LTG/2,:) = wall;

% diagonalizing the single-particle Hamiltonian
[ES,EV] = eig(TG_kin + diag(TG_trap));
TG_energies_ini = diag(EV);
TG_states = (1/sqrt(dx)).*ES;

% since T = 0 the chemical potential can be set freely 
Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1));
fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; 
fermi_ini = fermi_ini(fermi_ini>1e-8); % only consider occupied states 

% initial TG density from single-particle states
rho_ini = sum(fermi_ini.*abs(TG_states(:,1:length(fermi_ini))).^2,2);

clear TG_kin TG_trap ES EV TG_states Ef_ini fermi_ini

%% Ramping up the intra-species interaction strength
start = tic;
for i = 1:length(gMIXarr)
	gMIX = gMIXarr(i); fprintf('gMIX = %3.2f:\t',gMIX);
	tic; [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini); toc
	% using the previous ground state as the new initial state for speed-up
    wfi = d.wf;
	rho_ini = d.rho;
	TG_energies_ini = d.TG_energies;
	params(i) = p;
	data(i) = d;
	clear p d
    
    gap(i) = data(i).gap(end);
    mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    a0 = 0.5*gMIX^2/gBEC;  
    a0arr(i) = a0;
    eps = 6*a0^2/(5*mubar); 
    apinarr(i) = a0*(sqrt(1+2*eps)-1)/eps;
end
toc(start)

%% Write data file for densities in Fig. 1 (a) for gMIX = 0, 1, and 1.5 
mu0 = 0.5*gBEC*NBEC/posmax;
header = ["x","psi1","psi2","psi3","rho1","rho2","rho3"];
dataout = [x (abs(data(1).wf).^2 - mu0/gBEC) (abs(data(101).wf).^2 - mu0/gBEC) (abs(data(151).wf).^2 - mu0/gBEC) data(1).rho data(101).rho data(151).rho];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_NTG7_repulsive_densities.dat'),'Delimiter','tab');

%% Write data file for spectrum in Fig. 1 (b) for gMIX = 0 and 2
gMIX = 2;
mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
header = ["n","gm0","gm2"];
dataout = [(1:100).' data(1).TG_energies(1:100) (data(201).TG_energies(1:100) - gMIX*mubar/gBEC)];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_NTG7_gMIX2_spectrum.dat'),'Delimiter','tab');

%% Write data file for the numerical and analytical value of the energy gap as a function of gMIX in Fig. 1 (c)
header = ["gm","numerical","variational"];
dataout = [gMIXarr(2:end).' gap(2:end).' abs(apinarr(2:end).^2/6 - 2*a0arr(2:end).*apinarr(2:end)/3).'];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_NTG7_gMIXramp_gap.dat'),'Delimiter','tab');

clear start i wfi rho_ini TG_energies_ini mubar a0 eps gMIX header dataout
save(fname);