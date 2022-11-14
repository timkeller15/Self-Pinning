function [p,d] = selfpinning_groundstate(Ngrid,posmax,NTG,LTG,NBEC,gBEC,gMIX,T,wfi,rho_ini,TG_energies_ini)
    %% Simulation parameters
    % they are chosen to ensure a reasonably fast and accurate convergence
    % under most conditions but may be adjusted of course 
    
    dt = 1e-3;  % time-step for the imaginary time evolution
    cutoff = 1e-6;  % convergence criterion for the energy
    wall = 1e8; % box potential for the TG gas
    becsteps = 50; % number of steps for the imaginary time evolution after each TG gas density adjustment
    samples = 5e4; % maximum number of steps 
    [x,dx,px,~] = fftdef(posmax,Ngrid); % define position grid 

    % set initial densities and Tonks-Girardeau energies 
    wf = wfi;
    rho = rho_ini;
    TG_energies = TG_energies_ini;
    
    %% Bose-Einstein condensate
    Ekin = exp(-0.5*dt*px.^2); % kinetic energy operator for the BEC's imaginary time evolution
    Vtrap = zeros(Ngrid,1); % free space / no potential and periodic boundary conditions for the BEC 

    % initial BEC energy
    Ki = 0.5*conj(wfi).*ifft(px.^2.*fft(wfi));
    Vi = (Vtrap + gMIX*rho).*abs(wfi).^2 + 0.5*gBEC.*abs(wfi).^4;
    Ei = real(sum(Ki + Vi))*dx;
    
    %% Tonks-Girardeau gas 
    % set operators for single-particle Hamiltonian of the TG gas    
    TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
    TG_trap = zeros(Ngrid,1);
    TG_trap(abs(x)>LTG/2,:) = wall;
    
    % adjust the TG gases chemical potential to make sure that the particle
    % number is conserved
    Ef_ini = 0.5*(TG_energies_ini(NTG) + TG_energies_ini(NTG+1)); % initial guess for the chemical potential
    fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; 
    fermi_ini = fermi_ini(fermi_ini>1e-8);
    if T>0
        while abs(sum(fermi_ini) - NTG) > 1e-4
            Ef_ini = Ef_ini - T*(sum(fermi_ini) - NTG)/sum((exp((TG_energies_ini(1:length(fermi_ini)) - Ef_ini)./T).').*fermi_ini.^2); % Newton's method
            fermi_ini = 1./(exp((TG_energies_ini - Ef_ini)./T) + 1).'; 
            fermi_ini = fermi_ini(fermi_ini>1e-8); % only consider occupied states
        end
    end
    
    %% Initial values for the observables 
    energy = zeros(1,samples);
    gap = zeros(1,samples);
    Natom = zeros(1,samples);
    Nfermi = zeros(1,samples);
    f0 = zeros(1,samples);
    time = zeros(1,samples);
    
    energy(1) = Ei + sum(fermi_ini.'.*TG_energies_ini(1:length(fermi_ini))); 
    Natom(1) = sum(abs(wf).^2)*dx;
    Nfermi(1) = sum(fermi_ini);
    f0(1) = fermi_ini(1); 
    gap(1) = TG_energies(NTG+1) - TG_energies(NTG);

    %% Alternating imaginary time evolution for the BEC and diagonalization of the single-particle Hamiltonian for the TG gas until convergence 
    for i=2:samples
            
        % evolve the BEC wave function in imaginary time with a symmetric
        % split-step method
        for j = 1:becsteps
            % Update potential with TG gas density 
            V = Vtrap + gMIX*rho + gBEC*abs(wf).^2;
            
            % Step 1 in position space
            wf = exp(-0.5*dt*V).*wf;

            % Step 2 in momentum space
            wf = ifft(Ekin.*fft(wf));

            % Step 3 in position space
            wf = exp(-0.5*dt*V).*wf;

            % Normalize
            wf = wf/sqrt(sum(conj(wf).*wf)*dx/NBEC); 
        end
        
        % diagonalize the single-particle TG Hamiltonian with the updated BEC density 
        [ES,EV] = eig(TG_kin + diag(TG_trap + gMIX*abs(wf).^2));
        TG_energies = diag(EV);
        TG_states = (1/sqrt(dx)).*ES;

        % adjust the TG gases chemical potential to make sure that the particle
        % number is conserved
        Ef = 0.5*(TG_energies(NTG) + TG_energies(NTG+1)); % initial guess for the chemical potential
        fermi = 1./(exp((TG_energies - Ef)./T) + 1).'; 
        fermi = fermi(fermi>1e-8);
        if T>0
            while abs(sum(fermi) - NTG) > 1e-4
                Ef = Ef - T*(sum(fermi) - NTG)/sum((exp((TG_energies(1:length(fermi)) - Ef)./T).').*fermi.^2); % Newton's method 
                fermi = 1./(exp((TG_energies - Ef)./T) + 1).'; 
                fermi = fermi(fermi>1e-8); % only consider occupied states
            end
        end
        
        % calculate the updated TG gas density 
        rho = sum(fermi.*abs(TG_states(:,1:length(fermi))).^2,2);
            
        time(i) = becsteps*dt*(i - 1);

        % calculate observables 
        K = 0.5*conj(wf).*ifft(px.^2.*fft(wf));
        V = Vtrap.*abs(wf).^2 + 0.5*gBEC.*abs(wf).^4; % not counting Vint twice (only in TG gases energy)

        energy(i) = real(sum(K + V))*dx + sum(fermi.'.*TG_energies(1:length(fermi)));
        Natom(i) = sum(abs(wf).^2)*dx;
        Nfermi(i) = sum(fermi);
        f0(i) = fermi(1); 
        gap(i) = TG_energies(NTG+1) - TG_energies(NTG);

        % check for convergence 
        if i > 2
            diff = abs(energy(i) - energy(i-1));
            if diff < cutoff 
                time = time(1:i);
                energy = energy(1:i);
                Natom = Natom(1:i);
                Nfermi = Nfermi(1:i);
                f0 = f0(1:i);
                gap = gap(1:i);
                % return observables if convergence has been reached 
                p = v2struct(samples,becsteps,dt,cutoff,wall,posmax,Ngrid,Ei,Ef,NTG,LTG,NBEC,gBEC,gMIX,T);
                d = v2struct(wf,rho,TG_energies,energy,gap,time,Natom,Nfermi,f0,fermi);
                return
            end
        end 
    end

    % return observables if convergence hasn't been reached after the
    % maximum number of steps 
    p = v2struct(samples,becsteps,dt,cutoff,wall,posmax,Ngrid,Ei,Ef,NTG,LTG,NBEC,gBEC,gMIX,T);
    d = v2struct(wf,rho,TG_energies,energy,gap,time,Natom,Nfermi,f0,fermi);
end