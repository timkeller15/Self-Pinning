%% Load data created for Fig. 1
cd ..
fname = '/data/selfpinning_zero_temperature_NTG7_gMIXramp.mat';
fname = fullfile(fileparts(pwd),fname);
load(fname);

%% Extract the first/lowest eigenvalue of a pinned Tonks-Girardeau atom for the fit. The numerical values for the pinned eigenenergies are not perfectly identical, but using another pinned state or e.g. the mean only leads to a slight quantitative change of the resulting curves. 
for i = 1:length(data)
    gMIX = gMIXarr(i);
    mubar = 0.5*gBEC*(NBEC + gMIX*NTG/gBEC)/posmax;
    spectrum(i) = data(i).TG_energies(1) - gMIX*mubar/gBEC;
end
clear i gMIX mubar

%% Exclude gm = 0 points for which the model doesn't work
gMIXarr = gMIXarr(2:end).';
a0arr = a0arr(2:end); 
apinarr = apinarr(2:end);
spectrum = spectrum(2:end);

%% Calculate the different fits 
fit1 = abs(-0.5*a0arr.^2./spectrum - 1).';   % a0
fit2 = abs(-0.5*apinarr.^2./spectrum - 1).'; % apin
fit3 = abs((apinarr.^2/6 - 2*a0arr.*apinarr/3)./spectrum - 1).'; % variational
fit4 = abs(-0.5*a0arr.*apinarr./spectrum - 1).'; % geometric

%% Write data file for 
header = ["gm","a0","apin","variational","geometric"];
dataout = [gMIXarr fit1 fit2 fit3 fit4];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/selfpinning_spectrum_fit.dat'),'Delimiter','tab');