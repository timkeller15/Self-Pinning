%% Setting parameters
Ngrid = 2048;   % number of position grid points
posmax = 25;    % range of position grid
NTG = 7;    % number of atoms in the Tonks-Girardeau gas
LTG = 40;   % length of the box potential for the TG gas
NBEC = 1e4; % number of atoms in the BEC
gBEC = 1;   % BEC intra-species interaction strength
gMIX = 3;  % inter-species interaction strenght
temperatures = 0:.01:5; % temperature range
[x,dx] = fftdef(posmax,Ngrid);  % defines position grid
fname = '/data/selfpinning_finite_temperature_NTG7_gMIX3_momentum.mat';
fname = fullfile(fileparts(pwd),fname);
wall = 1e8;
[x,dx] = fftdef(posmax,Ngrid);

% set inital densities 
TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
TG_trap = zeros(Ngrid,1);
TG_trap(abs(x)>LTG/2,:) = wall;

[ES,EV] = eig(TG_kin + diag(TG_trap));
TG_energies_ini = diag(EV);
TG_states = (1/sqrt(dx)).*ES;

T = 0;
Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1));
fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; 
fermi_ini = fermi_ini(fermi_ini>1e-8);
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
    
	% find the ground-state first
    [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini);
	wfi = d.wf;
	rho_ini = d.rho;
	TG_energies_ini = d.TG_energies;
	params(i) = p;
	data(i) = d;
	clear p d
    
    % calculate the rspdm 
    tic; rspdm = selfpinning_rspdm(Ngrid,posmax,LTG,gMIX,temperatures(i),params(i).Ef,data(i).wf); toc
    
    [ES,EV] = eig(rspdm);
    occupationnr = flipud(diag(EV))*dx;
    naturalorbitals=(1/sqrt(dx))*fliplr(ES);
    occupationnr = occupationnr(occupationnr>0);
    naturalorbitals = naturalorbitals(:,occupationnr>0);
    
    % calculate the momentum distribution from the diagonalized rspdm 
    momdis(:,i) = fftshift(abs(fft(naturalorbitals)*(dx/sqrt(2*pi))).^2*occupationnr); % normalisation factor / Parseval's theorem 
    clear ES EV rspdm occupationnr naturalorbitals
    save(fname);
end
clear i wfi rho_ini TG_energies_ini T  
save(fname);

%% Calculating the model data for the momentum distribution at gMIX = 3
% For better alignment with the numerical data we use the numerically
% obtained spectrum for the continuum states needed to determine the
% chemical potential. However, a complete ab-initio approach with the
% continuum spectrum given by the unperturbed box potential according to
% continuum = 0.5*(pi/LTG)^2*(1:(Ngrid-NTG)).^2; 
% also yields good agreement. 

gMIX = 3; 
mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
continuum = (data(1).TG_energies(NTG+1:end) - gMIX*mubar/gBEC).';

a0 = 0.5*gMIX^2/gBEC;
for i = 1:length(temperatures)
    % reset inital guesses 
    T = temperatures(i); 
    fpin = 1; 
    fpin_old = 0; 
    count = 0;
    mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    eps = 6*a0^2/(5*mubar); 
    
    % solve the self-consistency equation for the Fermi-Dirac factors 
    while abs(fpin-fpin_old) > 1e-4 
        fpin_old = fpin; 
        Epin = -0.5*a0^2*(sqrt(1 + 2*eps*fpin^2) - 1)/eps;
        TG_energies_ana = [Epin*ones(1,NTG) continuum];
        Ef = 0.5*(TG_energies_ana(NTG+1) + TG_energies_ana(NTG));
        fermi = 1./(exp((TG_energies_ana - Ef)./T) + 1); 
        fermi = fermi(fermi>1e-8);
        iter = 0;
        % adjust the chemical potential to conserve the particle number
        while abs(sum(fermi) - NTG) > 1e-4
            iter = iter + 1;
            Ef = Ef - T*(sum(fermi) - NTG)/sum(exp((TG_energies_ana(1:length(fermi)) - Ef)./T).*fermi.^2);
            fermi = 1./(exp((TG_energies_ana - Ef)./T) + 1); fermi = fermi(fermi>1e-8);
             if iter > 1e3
                fpin = NaN; 
                break
            end
        end
        
        fpin = fermi(1);  % update the Fermi-Dirac factors 
        mubar = 0.5*gBEC*(NBEC + gMIX*sqrt(fpin)*NTG/gBEC)/posmax;
        eps = 6*a0^2/(5*mubar); 
        count = count + 1; 
        if count > 1e3
            fpin = NaN; 
            break
        end
    end
    pinnedness_ana(i) = fpin; 
    apinarr(i) = fpin*a0*(sqrt(1 + 2*eps*fpin^2) - 1)/(eps*fpin^2);
end
clear i a0 T fpin fpin_old count mubar eps Epin TG_energies_ana Ef fermi iter 

% zero-momentum peaks according to the analytical model 
mompeak_ana = 0.25*pi./apinarr;

a0 = 0.5*gMIX^2/gBEC; mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax; eps = 6*a0^2/(5*mubar); ap = a0*(sqrt(1+2*eps)-1)/eps;
momdis_pin = (0.25*pi/ap)*1./cosh(0.5*pi*k/ap).^2; 

ind = find(islocalmax(momdis(Ngrid/2 + 1,:)),1,'last'); % find the position of the last 'de-pinning' / first data point in the thermal state. might require manual adjusting
sigma = 1/(sqrt(2*pi)*momdis(Ngrid/2 + 1,ind));
momdis_gauss = exp(-0.5*(k/sigma).^2)/(sigma*sqrt(2*pi));

%% Write data file for the momentum distributions in Fig. 3 (a) and (b)
header = ["k", "momdis_zero", "momdis_pin",	"momdis_thermal", "momdis_gauss"];
dataout = [k.' momdis(:,1) momdis_pin.' momdis(:,ind) momdis_gauss.'];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_NTG7_gMIX3_momentum_distribution.dat'),'Delimiter','tab');

%% Write data file for the zero-momentum peaks in Fig. 3 (c)
Tf = 0.5*a0*ap; % scale temperature in units of the energy gap size at T = 0 
header = ["T","mompeak_num","mompeak_ana"];
dataout = [temperatures.'/Tf momdis(Ngrid/2 + 1,:).' mompeak_ana.'];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_NTG7_gMIX3_momentum_peaks.dat'),'Delimiter','tab');

clear a0 mubar eps ap header dataout
save(fname);