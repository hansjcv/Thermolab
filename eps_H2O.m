function eps_di = eps_H2O(T,P,rho,dielectric_model)
% Preprocess
P = P/1e8;
switch dielectric_model
    case 'JN91'
        e =   [
            0.1470333593e2,  0.2128462733e3, -0.1154445173e3, 0.1955210915e2, -0.8330347980e2, ...
            0.3213240048e2, -0.6694098645e1, -0.3786202045e2, 0.6887359646e2, -0.2729401652e2
            ];
        rho_gcm3 = rho/1e3;
        Tref     = 298.15;
        Tn       = T/Tref;
        kT_1     = ones(size(Tn));
        kT_2     = e(1)./Tn;
        kT_3     = e(2)./Tn    + e(3)     + e(4)*Tn;
        kT_4     = e(5)./Tn    + e(6)*Tn  + e(7)*Tn.^2;
        kT_5     = e(8)./Tn.^2 + e(9)./Tn + e(10);
        eps_di   = kT_1.*rho_gcm3.^0 + kT_2.*rho_gcm3.^1 + kT_3.*rho_gcm3.^2 + kT_4.*rho_gcm3.^3 + kT_5.*rho_gcm3.^4;
    case 'F97'
        T        = T;
        rho_kgm3 = rho;
        % Constants
        eps_0 = (4e-7*pi*(299792458)^2)^(-1);
        k_B   = 1.380658e-23;
        Na    = 6.0221367e23;
        Mw    = 0.018015268;
        alph  = 1.636e-40;
        mu    = 6.138e-30;
        Tc    = 647.096;
        % Preprocess
        rho_c = 322/Mw;% mol/m3
        rho   = rho_kgm3/Mw;%mol/m3
        % Parameters
        N = [
            0.978224486826
            -0.957771379375
            0.237511794148
            0.714692224396
            -0.298217036956
            -0.108863472196
            0.949327488264e-1
            -0.980469816509e-2
            0.16516763497e-4
            0.937359795772e-4
            -0.12317921872e-9
            0.196096504426e-2];
        i = [1 1 1 2 3 3 4 5 6 7 10];
        j = [0.25 1 2.5 1.5 1.5 2.5 2 2 5 0.5 10];
        q = 1.2;
        % Process
        g = 1;
        for k = 1:11
            g = g + N(k)*(rho/rho_c).^(i(k)).*(Tc./T).^(j(k));
        end
        g       = g + N(12)*(rho/rho_c).*(T/228 - 1).^(-q);
        A       = (Na*mu^2*rho.*g)./(eps_0*k_B*T);
        B       = (Na*alph*rho)./(3*eps_0);
        eps_di = (1 + A + 5*B + sqrt(9 + 2*A + 18*B + A.^2 + 10*A.*B + 9*B.^2))./(4 - 4*B);
    case 'S14'
        rho_gcm3 = rho/1e3;
        TC = T - 273.15;
        %Relevant parameters
        a1 = -1.57637700752506e-03;
        a2 = 6.81028783422197e-02;
        a3 = 0.754875480393944;
        b1 = -8.01665106535394e-05;
        b2 = -6.87161761831994e-02;
        b3 = 4.74797272182151;
        
        A = a1 * TC + a2 * sqrt(TC) + a3;
        B = b1 * TC + b2 * sqrt(TC) + b3;
        
        eps_di = exp(B) .* rho_gcm3.^A; % equation 8 in Sverjensky GCA 2014
end
end