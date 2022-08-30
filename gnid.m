function g_nid = gnid(T,P,p,model_type,v0,alp,w,rho_w,eps_di,chg)
P = P/1e8;
if model_type == 1 || model_type == 4
    alp  = [1 T P]*alp;
    W    = w(:,:,1)   +   w(:,:,2)*T +   w(:,:,3)*P;
    g_nid = 0;
    for i = 1:size(W,1)
        for j = 1:size(W,2)
            g_nid = g_nid + p(:,i).*p(:,j).*alp(i).*alp(j)./(p*alp')./(alp(i)+alp(j)).*W(i,j);
        end
    end
elseif model_type == 2    
    z          = p'; % for this model type, the site fractions are equal to proportions.
    rho_w_gcm3 = rho_w/1e3;
    % HKF
    R     = 8.3145;
    a0    = 3.72e-8; % Angstrom a0 parameter, fixed now
    b_gam = 0.03; % fixed b gamma parameter also used in neutral species
    nu    = 2; % stoichiometric coefficient of charged species in the fluid
    Mw    = 18.01528e-3;
    Nw    = 1/Mw;
    l_gam = zeros(size(z,1),size(z,2));
    I  = 0;
    mz = 0;
    for i = 1:size(z,1)
        z_m(i,:)     = z(i,:)./(z(end,:)*Mw); % Try again
    end
    for i = 1:size(z,1) % loop over species in the fluid
        I      = I + 0.5*z_m(i,:)*chg(i).^2; % True Ionic strength
        if chg(i)~=0
            mz = mz + z_m(i,:); % total molality of charged species
        end
    end
    mt   = sum(z_m(1:end-1,:),1); % total molality
    Lgam = -log10(1+Mw*mt); % large gamma (for 1 kg of water) conversion factor see Helgeson L gamma page 1293 eq 122
    % Debye-Hueckel A and B parameters
    A = (1.82483e6)*sqrt(rho_w_gcm3(:)) ./ (T(:).*eps_di(:)).^(3/2);
    B = (50.2916e8)*sqrt(rho_w_gcm3(:)) ./ sqrt(T(:).*eps_di(:));
    for iPT = 1:length(T(:))
        % log Activity coefficient of solutes
        for i = 1:size(z,1)-1
            l_gam(i,:) = - (A(iPT) * chg(i)^2 * sqrt(I))./(1+a0*B(iPT)*sqrt(I)) + b_gam*I + Lgam; % probably molar scale see Helgeson L gamma page 1293 eq 122
        end
        % Activity of water
        lam       = 1 + a0*B(iPT)*sqrt(I); % lamda
        sig       = 3./(a0^3*B(iPT).^3*sqrt(I.^3)).*(lam-1./lam-2*log(lam));  % sigma coefficient
        phi       = -log(10)*mz./mt.*(A(iPT).*sqrt(I).*sig/3 + Lgam./(Mw*nu*I) - b_gam*I/2); % Osmotic coefficient
        phi(I==0) = 0; % set osmotic coefficient to 0 for case of no electrolytes
        a_w       = exp(-phi.*mt/Nw); % Activity of water
        gam_w     = a_w./z_m(end,:); % Activity coefficient water DIVIDED by z_molality 26-4-2022       
%         gam_w     = a_w./z(end,:); % Activity coefficient water DIVIDED by z mole fraction
        % Activity coefficients
        gam     = 10.^(l_gam);
        % Make the excess g in mole fractions
        g_ex   = zeros(1,size(l_gam,2));
        for i = 1:size(z,1)-1
            g_ex   = g_ex  + R*T(iPT)*z(i,:).*(log(gam(i,:)) + log(Nw)); % mole fraction interpretation
        end
        % and add the water contribution
        g_ex = g_ex + R*T(iPT)*z(end,:).*log(gam_w);
    end
    g_nid = g_ex'; 
elseif model_type == 3
    V = v0;
    % Mixing parameters1
    w_data   = [12893, -6.501,1.0112];
    A = w_data(1);
    B = w_data(2);
    C = w_data(3);
    W_AN99 = (A+B*T).*(1-exp(-20*P)) + C*T.*P;
    w = [0 W_AN99/2
        W_AN99/2 0];
    g_nid = 0;
    for i = 1:size(w,1)
        for j = 1:size(w,2)
            g_nid = g_nid + p(:,i).*p(:,j).*V(i).*V(j)./(p*V')./(V(i)+V(j)).*w(i,j);
        end
    end
elseif model_type == 5 % Evans and Powell (2006) unfinished model!
      g_nid = sum(p(:,1:end-1)*8.3145*T*log(1000/18.0150),2);
elseif model_type == 0
    g_nid = 0;
end
end