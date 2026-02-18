function rho_w = rho_H2O(T,P,EOS)
P     = P/1e8;
rho_w = zeros(length(T),1);
switch EOS
    case 'JN91'
        %% Johnson and Norton 1991
        rho_g1 = 500;
        rho_g2 = 1000;
        rho_g3 = 1600;
        for i = 1:length(T)            
            if P(i)<2
                if T(i) > 573.15
                    rho_w(i) = fzero(@(rho) JN91(rho,T(i))-P(i),rho_g1); % T in K , P in kb
                else
                    rho_w(i) = fzero(@(rho) JN91(rho,T(i))-P(i),rho_g2); % T in K , P in kb
                end
            else
                rho_w(i) = fzero(@(rho) JN91(rho,T(i))-P(i),rho_g3); % T in K , P in kb
            end
        end
    case 'IAPWS'
        %% IAPWS EOS
        P  = P*1e2;
        R  = 0.46151805e3; Tc = 647.096; Pc = 22.064; rhoc = 322;
        a1 = -7.85951783;  a2 = 1.84408259; a3 = -11.7866497;a4 = 22.6807411;a5 = -15.9618719; a6 = 1.80122502;
        b1 = 1.99274064;   b2 = 1.09965342; b3 = -0.510839303; b4 = -1.75493479; b5 = -45.5170352; b6 = -6.74694450e5;
        c1 = -2.03150240; c2 = -2.68302940; c3 = -5.38626492; c4 = -17.2991605;c5 = -44.7586581; c6 = -63.9201063;
        max_err = 1e-13;
        for iT = 1:length(T)
            thet  = 1-T(iT)/Tc;
            Ps    = Pc*exp(Tc/T(iT)*(a1*thet + a2*thet^(1.5) + a3*thet^3 + a4*thet^(3.5) + a5*thet^4 + a6*thet^(7.5)));
            rhof  = rhoc*(1 + b1*thet^(1/3) + b2*thet^(2/3) + b3*thet^(5/3) + b4*thet^(16/3) + b5*thet^(43/3) + b6*thet^(110/3));
            rhov  = rhoc*exp( c1*thet^(2/6) + c2*thet^(4/6) +c3*thet^(8/6)+c4*thet^(18/6)+c5*thet^(37/6) + c6*thet^(71/6));
            p1100 = (24.60396371*1e5*T(iT) -4695.30045565*1e5)/1e6;
            p800  = (13.7263594135030e5*T(iT) -7130.8669147992003e5)/1e6;
            if P(iT)>p1100
                rho_g1 = 1200;
            elseif P(iT)<Pc
                if P(iT)<=0.9999*Ps
                    if P(iT) < 0.05
                        rho_g1 = 0.1;
                    elseif P(iT) > 10
                        rho_g1 = 100;
                    else
                        rho_g1 = 1;
                    end
                elseif P(iT)>=1.00005*Ps
                    if P(iT) < p800
                        rho_g1 = 800;
                    else
                        rho_g1 = 1000;
                    end
                else
                    Ps0 = Ps;
                    delP = 1;
                    if T(iT) <= Tc - 3.5e-3 %If T is not close to Tc,
                        while delP > max_err
                            [~,~,~,~,phir_f] = IAPWS(rhof,T(iT));
                            [~,~,~,~,phir_v] = IAPWS(rhov,T(iT));
                            Ps   = R*T(iT)*((phir_f-phir_v) + log(rhof/rhov)) / (1/rhov - 1/rhof);
                            rhof = fzero(@(rho) IAPWS(rho,T(iT))-Ps/1e6,rhof);
                            rhov = fzero(@(rho) IAPWS(rho,T(iT))-Ps/1e6,rhov);
                            delP = abs((Ps - Ps0)/Ps);
                            Ps0 = Ps;
                        end
                    end
                    if T(iT) <= Tc - 6e-6
                        rhoMin = fminbnd(@(rho)  IAPWS(rho,T(iT)), rhov, rhof);
                        rhoMax = fminbnd(@(rho) -IAPWS(rho,T(iT)), rhov, rhof);
                        drho = 1;
                        while delP > max_err
                            [~,~,~,~,phir_f] = IAPWS(rhof,T(iT));
                            [~,~,~,~,phir_v] = IAPWS(rhov,T(iT));
                            Ps   = R*T(iT)*((phir_f-phir_v) + log(rhof/rhov)) / (1/rhov - 1/rhof);
                            rhof = fzero(@(rho) IAPWS(rho,T(iT))-Ps/1e6,[rhoMin, rhof + drho]);
                            rhov = fzero(@(rho) IAPWS(rho,T(iT))-Ps/1e6,[rhov - drho, rhoMax]);
                            delP = abs((Ps - Ps0)/Ps);
                            Ps0 = Ps;
                        end
                    end
                    if P <= Ps/(1 + 5*eps)    % steam
                        rho_g1 = rhov;
                    elseif P >= Ps*(1 + 5*eps)   % water
                        rho_g1 = rhof;
                    end
                end
            else
                rho_g1 = 1000;
            end
            if P(iT)>90
                rho_g1 = 1500;
            end
            rho_w(iT)  = fzero(@(rho) IAPWS(rho,T(iT))-P(iT),rho_g1);
        end        
    case 'ZD09'
        %% Zhang and Duan 2009 EOS
        rho_g1 = 1800;
        rho_g2 = 1000;
        for i = 1:length(T)
            if T(i) > 673.15
                rho_w(i) = fzero(@(rho) rho_ZD09(rho,T(i),P(i)),rho_g1); % T in K , P in kbar
            else
                warning('Out of calibration range for Zhang and Duan EOS, switching to IAPWS');
                if P(i)<2
                    rho_w(i) = fzero(@(rho) rho_IAPWS(rho,T(i),P(i)*1e2),rho_g2); % T in K , P in MPa
                else
                    rho_w(i) = fzero(@(rho) rho_IAPWS(rho,T(i),P(i)*1e2),rho_g1); % T in K , P in MPa
                end
            end
            
        end
    case 'ZD05'
        %% Zhang and Duan 2005 EOS
        rho_g1 = 1800;
        rho_g2 = 1000;
        for i = 1:length(T)            
            if T(i) > 573.15
                rho_w(i) = fzero(@(rho) ZD05(rho,T(i))-P(i),rho_g1); % T in K , P in kbar
            else
                rho_w(i) = fzero(@(rho) ZD05(rho,T(i))-P(i),rho_g2); % T in K , P in kbar
            end
        end
    case 'CORK'
        %% Holland and Powell (1991) CORK
        td = [-1.722516050000000, 0.001888000000000, 0,0.000401000000000,   0.000000086560000,   4.875000000000000,  -0.002512000000000,                   0,                   0,                   0,                   0,                   0,                   0,                   0,                   0,                   0,   0.010000000000000,                   0]*1e2;
        [~,V] = intVdP(T,P*1e8,td,4);
        rho_w = (18.0150e3/10)./V;
    case 'PS94'
        td = [-1.722516050000000,0.001888000000000, 0,0.000401000000000,0.000000086560000,4.875000000000000, -0.002512000000000,0, 0,0, 0,NaN,NaN, NaN,NaN,NaN, NaN, NaN, NaN, NaN, NaN, NaN,0.030000000000000,  0, 0.010000000000000,   0, 0]*1e2;
        [~,V] = intVdP(T,P*1e8,td,5);
        rho_w = (18.0150e3/10)./V;
end
end