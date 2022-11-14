%% Setting parameters
Ngrid = 2048;   % number of position grid points
posmax = 25;    % range of position grid
NTG = 7;    % number of atoms in the Tonks-Girardeau gas
LTG = 40;   % length of the box potential for the TG gas
NBEC = 1e4; % number of atoms in the BEC
gBEC = 1;   % BEC intra-species interaction strength
[x,dx] = fftdef(posmax,Ngrid);  % defines position grid
wall = 1e8; % box potential height
temperatures = 0:.01:5; % temperature range
gMIXarr = -5:.01:5; % interspecies interaction range

%% Calculating the ground-states. 
% The outer loop over the interaction strenghts should be parallelized as
% much as possible for reasonable execution times. For example by commenting it out
% and setting gMIX = gMIXarr(${SLURM_ARRAY_TASK_ID}) when running the code
% as an array job on a cluster using slurm. 
% It is also sufficient to only consider gMIXarr = 0:.01:5 since we have
% checked that the resulting states for the Tonks-Girardeau gas are
% symmetric with respect to the sign of gMIX. 

pinnedness_num = zeros(length(gMIXarr),length(temperatures));
start = tic;
for i = 1:length(gMIXarr)
    gMIX = gMIXarr(i); 
    fprintf('gMIX = %3.2f:\t',gMIX);
    fname = sprintf('/data/selfpinning_finite_temperature_ramp_NTG7_gMIX%3.4f.mat',gMIX); 
    fname = fullfile(fileparts(pwd),fname);
    
    % resetting the initial conditions for each iteration of gMIX
    T = 0; 
    
    TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
    TG_trap = zeros(Ngrid,1);
    TG_trap(abs(x)>LTG/2,:) = wall;

    [ES,EV] = eig(TG_kin + diag(TG_trap));
    TG_energies_ini = diag(EV);
    TG_states = (1/sqrt(dx)).*ES;

    Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1));
    fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; 
    fermi_ini = fermi_ini(fermi_ini>1e-8);
    rho_ini = sum(fermi_ini.*abs(TG_states(:,1:length(fermi_ini))).^2,2);

    wfi = sqrt(NBEC/(Ngrid*dx))*ones(Ngrid,1);
    clear TG_kin TG_trap ES EV TG_states fermi_ini Ef_ini

    % ramping up the temperature 
    tic
    for j = 1:length(temperatures)
        T = temperatures(j); 
        [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini);
        wfi = d.wf;
        rho_ini = d.rho;
        TG_energies_ini = d.TG_energies;
        params(j) = p;
        data(j) = d;
        clear p d
        
        pinnedness_num(i,j) = sum(data(j).fermi(1:NTG))/NTG;
    end
    toc
    clear j T wfi rho_ini TG_energies_ini
    save(fname, '-regexp', '^(?!pinnedness_num$).');
    clear params data
end
toc(start)
clear start i gMIX fname 

%% Calculating the model data for the ground state occupancy at gMIX = 3
% For better alignment with the numerical data we use the numerically
% obtained spectrum for the continuum states needed to determine the
% chemical potential. However, a complete ab-initio approach with the
% continuum spectrum given by the unperturbed box potential according to
% continuum = 0.5*(pi/LTG)^2*(1:(Ngrid-NTG)).^2; 
% also yields good agreement. 

gMIX = 3; 
mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
fname = sprintf('/data/selfpinning_finite_temperature_ramp_NTG7_gMIX%3.4f.mat',gMIX); 
fname = fullfile(fileparts(pwd),fname);
load(fname,'data');
continuum = (data(1).TG_energies(NTG+1:end) - gMIX*mubar/gBEC).';
clear fname data

a0 = 0.5*gMIX^2/gBEC;
for i = 1:length(temperatures)
    T = temperatures(i); 
    fpin = 1; 
    fpin_old = 0; 
    count = 0;
    mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    eps = 6*a0^2/(5*mubar); 
    
    while abs(fpin-fpin_old) > 1e-4 
        fpin_old = fpin; 
        Epin = -0.5*a0^2*(sqrt(1 + 2*eps*fpin^2) - 1)/eps;
        TG_energies_ana = [Epin*ones(1,NTG) continuum];
        Ef = 0.5*(TG_energies_ana(NTG+1) + TG_energies_ana(NTG));
        fermi = 1./(exp((TG_energies_ana - Ef)./T) + 1); 
        fermi = fermi(fermi>1e-8);
        iter = 0;
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
end
clear i gMIX a0 T fpin fpin_old count mubar eps ap TG_energies_ana Ef fermi iter 


%% Calculating the analytical line for the critical temperature
for i = 1:length(gMIXarr)
    gMIX = gMIXarr(i); 
    a0 = 0.5*gMIX^2/gBEC;
    mubar0 = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    eps = 6*a0^2/(5*mubar0); 
    syms u; fstars = double(vpasolve(4*eps*u^3 - 2*eps*u^2 + 3*u - 2 == 0,u,2/3)); 
    fstar = fstars(1); clear u; clear fstars
    E = -0.5*a0^2*(sqrt(1+2*eps*fstar^2) - 1)/eps;
    Tcrit(i) = 3*E/pi^2*log(1/fstar-1)*(sqrt(1 + pi^2/(3*log(1/fstar-1)^2))-1);
end
clear i gMIX a0 mubar0 eps E fstar

fname = '/data/selfpinning_finite_temperature_phasediagram.mat';
fname = fullfile(fileparts(pwd),fname);
save(fname);

%% Write data file for the numerical and analytical value of the ground state occupancy as a function of temperature in Fig. 2 (a)
gMIX = 3; a0 = 0.5*gMIX^2/gBEC; mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax; eps = 6*a0^2/(5*mubar); ap = a0*(sqrt(1+2*eps)-1)/eps;
Tf = 0.5*a0*ap; % scale temperature in units of the energy gap size at T = 0
header = ["T","PN_num","PN_ana"];
dataout = [temperatures.'/Tf pinnedness_num(gMIXarr==3,:).' pinnedness_ana.']; 
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_NTG7_gMIX3_pinnedness.dat'),'Delimiter','tab');

%% Create heatmap plot for Fig. 2 (b)
load(fullfile(fileparts(pwd),'/data/smoothcoolwarm.mat'));
figure
set(gca,'position',[0 0 1 1]); hold all
surf(temperatures,gMIXarr,pinnedness_num); 
shading flat; view(2);
colormap(flipud(smoothcoolwarm));
caxis([0 1])
saveas(gcf,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_NTG7_phasediag.png'))

%% Write data file for the critical temperature line in Fig. 2 (b)
header = ["Tcrit","gMIX"];
dataout = [Tcrit.'/Tf gMIXarr.']; % scale temperature in units of the energy gap size at T = 0 for fixed gMIX = 3
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_finite_temperature_NTG7_Tcrit.dat'),'Delimiter','tab');
