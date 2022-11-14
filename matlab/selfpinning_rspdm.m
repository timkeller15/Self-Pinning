function rspdm = selfpinning_rspdm(Ngrid,posmax,LTG,gMIX,T,Ef,wf)
    % calculate single-particle states for the immersed Tonks-Girardeau gas in a box
    % potential for given BEC wave function wf.
    [x,dx] = fftdef(posmax,Ngrid);
    TG_kin = -1/(2*dx^2)*(diag(ones(Ngrid-1,1),-1) - 2*diag(ones(Ngrid,1)) + diag(ones(Ngrid-1,1),1));
    TG_trap = zeros(Ngrid,1); 
    wall = 1e8; 
    TG_trap(abs(x)>LTG/2,:) = wall;
    [ES,EV] = eig(TG_kin + diag(TG_trap + gMIX*abs(wf).^2));
    TG_energies = diag(EV); 
    TG_states = (1/sqrt(dx)).*ES; 
    clear ES EV
    
    % calculate the reduced single-particle density matrix for the TG gas
    % according to Atas et al., Phys. Rev. A 95, 043622 (2017).
    rspdm = zeros(Ngrid);
    fermi = 1./(exp((TG_energies - Ef)./T) + 1).';
    % only consider states with finite occupation 
    fermi = fermi(fermi>1e-8);
    phiP = sqrt(2*dx*fermi).*TG_states(:,1:length(fermi));
    phiQ = sqrt(fermi).*TG_states(:,1:length(fermi));
    for i = 1:Ngrid 
        bigP = eye(length(fermi));
        for j = i:-1:1 % only calcuate lower triangle since the rspdm is symmetric 
            bigP = bigP - (phiP(j,:)'*phiP(j,:)); 
            if i==j
                rspdm(i,j) = phiQ(j,:)*phiQ(i,:)';
            else
                bigQ = inv(bigP).'.*det(bigP);
                rspdm(i,j) = phiQ(j,:)*bigQ*phiQ(i,:)'; 
            end
        end
    end
    % mirror the full rspdm from lower triangle 
    rspdm = (tril(rspdm,-1) + tril(rspdm,-1).' + diag([0 diag(rspdm,-1).']))/sum(fermi);
end