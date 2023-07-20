function [gam,a_w] = gam_HKF(T,z,chg,rho_w,eps_di)
rho_w_gcm3 = rho_w/1e3;
% HKF
R     = 8.3145;
a0    = 4.5e-8; % a0    = 3.72e-8; % Angstrom a0 parameter, fixed now
b_gam = 0.03; % fixed b gamma parameter also used in neutral species
nu    = 2; % stoichiometric coefficient of charged species in the fluid
Mw    = 18.01528e-3;
Nw    = 1/Mw;
l_gam = zeros(size(z,1),size(z,2));
I  = 0;
mz = 0;
z_m = z;
for i = 1:size(z,1) % loop over species in the fluid
    I      = I + 0.5*z_m(i,:)*chg(i).^2; % True Ionic strength
    if chg(i)~=0
        mz = mz + z_m(i,:); % total molality of charged species
    end
end
mt   = sum(z_m,1); % total molality
Lgam = -log10(1+Mw*mt); % large gamma (for 1 kg of water) conversion factor see Helgeson L gamma page 1293 eq 122
% Debye-Hueckel A and B parameters
A = (1.82483e6)*sqrt(rho_w_gcm3(:)) ./ (T(:).*eps_di(:)).^(3/2);
B = (50.2918649e8)*sqrt(rho_w_gcm3(:)) ./ sqrt(T(:).*eps_di(:));
for iPT = 1:length(T(:))
    % log Activity coefficient of solutes
    for i = 1:size(z,1)
        l_gam(i,:) = - (A(iPT) * chg(i)^2 * sqrt(I))./(1+a0*B(iPT)*sqrt(I)) + b_gam*I + Lgam; % probably molar scale see Helgeson L gamma page 1293 eq 122
    end
    % Activity of water
    lam       = 1 + a0*B(iPT)*sqrt(I); % lamda
    sig       = 3./(a0^3*B(iPT).^3*sqrt(I.^3)).*(lam-1./lam-2*log(lam));  % sigma coefficient
    phi       = -log(10)*mz./mt.*(A(iPT).*sqrt(I).*sig/3 + Lgam./(Mw*nu*I) - b_gam*I/2); % Osmotic coefficient
    phi(I==0) = 0; % set osmotic coefficient to 0 for case of no electrolytes
    a_w       = exp(-phi.*mt/Nw); % Activity of water
    % Activity coefficients
    gam     = 10.^(l_gam);
end