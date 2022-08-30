function [VdP,V] = intVdP(T,P,td,eos)
P  = P/1e8;
if eos==4
    fl = td(17);
    [VdP,V] = HP_CORK(T,P,fl);
    VdP = VdP*1e3;
elseif eos == 5
    fl = td(25);
    [VdP,V] = HP_CORK(T,P,fl);
    VdP = VdP*1e3;
elseif eos == 1
    Tr = 25 + 273.15;
    Pr = 1e-3;Vr = td(3);a0 = td(8);k0 = td(10);dqf = td(18);
    V0  = Vr*(1  + a0*(T-Tr)-20*a0*(sqrt(T)-sqrt(Tr)));% Volume at constant P
    kT  = k0*(1-1.5e-4*(T-Tr)); % Compressibility at constant P
    V   = V0.*(1 + 4*(P-Pr)./kT).^(-1/4); % Volume
    VdP = V0.*kT/3.*((1+4.*(P-Pr)./kT).^(3/4)-1);% Integral of volume
    if dqf ==1
        VdP = Vr*P;
    end
    VdP = VdP*1e3;
elseif eos == 2
    Tr = 25 + 273.15;
    Pr = 1e-3;
    %% Load data
    Sref = td(2);V1_298 = td(3); a0   = td(8);
    k0     = td(9);k0p    = td(10);k0pp = td(11);
    dk0dt = td(22); sum_nphs = td(23); isL = td(24);dqf = td(27);
    %% TEOS (see Holland and Powell 2011)
    theta = 10636/((Sref*1e3)/sum_nphs + 6.44);
    u     = theta./T;
    u0    = theta./Tr;
    ksi0  = (u0.^2.*exp(u0))./((exp(u0)-1).^2);
    if isL == 1
        pth   = 0;
        kT    = k0 + dk0dt.*(T-Tr);
        V1_T  = V1_298.*exp(a0.*(T-Tr));
    else
        pth   = a0.*k0.*(theta./ksi0).*(1./(exp(u)-1)-1./(exp(u0)-1));
        kT    = k0;
        V1_T  = V1_298;
    end
    ta    = (1+k0p)./(1+k0p+kT*k0pp);
    tb    = k0p./kT - k0pp/(1 + k0p);
    tc    = (1 + k0p + kT .* k0pp)./(k0p + k0p.^2 - kT .* k0pp);
    VdP = V1_T .* ((P - Pr).*(1 - ta) + ta.*(-(1 + tb .* (P - pth)).^(1 - tc) ...
        +           (1 + tb .* (Pr - pth)).^(1 - tc))./(tb.*(tc - 1)))...
        ./         ((1 - ta) + ta.*(1 + tb * Pr).^(-tc));
    VdP(isnan(VdP)) = 0;
    if dqf == 1
        VdP = V1_T.*P;
    end
    VdP = VdP*1e3;
    V   = V1_T;
elseif eos == 3 % HKF
    Pref   = 1;
    Pbar   = P*1e3;
    psi    = 2.6e3;theta = 228;td(4)  = td(4)*1e-1;td(5)  = td(5)*1e2;td(7) = td(7)*1e4; td(9)  = td(9)*1e4;td(10) = td(10)*1e5;
    % HKF standard molal gibbs energies
    V   = td(4) + td(5)./(psi + Pbar) + td(6)./(T - theta) + td(7)./((psi + Pbar).*(T - theta));
    VdP = (td(4)*(Pbar-Pref) ...
        +  td(5)*log((psi+Pbar)/(psi+Pref)) ...
        + (1./(T-theta)) .* (td(6)*(Pbar-Pref) ...
        +                    td(7)*log((psi+Pbar)./(psi+Pref))))*4.184; %derived from eq. 59 Johnson et al 1992
elseif eos == 10 % Helgeson, Supcrt Minerals
    Tr = 298.15;
    Pr = 1;
    P  = P*1e3;
    v0r = td(4);a = td(5);b = td(6);c = td(7);
    Ttr1     = td(8); Htr1  = td(9);Vtr1 = td(10); dPdTtr = td(11); atr = td(12); btr = td(13); ctr = td(14);
    Ttr2    = td(15);Htr2 = td(16);Vtr2 = td(17);dPdTtr2 = td(18);atr2 = td(19);btr2 = td(20);ctr2 = td(21);
    Ttr3    = td(22);Htr3 = td(23);Vtr3 = td(24);dPdTtr3 = td(25);atr3 = td(26);btr3 = td(27);ctr3    = td(28);
    tr_id   = td(30);  itype   = td(31);
    Ttr = [Tr Ttr1 Ttr2 Ttr3 ];Htr = [0  Htr1 Htr2 Htr3 ];Vtr = [v0r Vtr1 Vtr2 Vtr3 ];atr = [a   atr atr2 atr3 ];btr = [b   btr btr2 btr3 ];ctr = [c   ctr ctr2 ctr3 ];
    VdP    = zeros(size(T'));
    V      = ones(size(T'))*v0r;
    for iPT = 1:length(T)
        Gpterm = 0;
        Tlo = Ttr(abs(atr)>0);
        for i = 1:tr_id+1
            T1 = Tlo(i);
            if T1>T(iPT), break,end
            VdP(iPT) = VdP(iPT) + Vtr(i)*(P(iPT) - Pr);
            V(iPT)   = V(iPT)   + Vtr(i);
            Gpterm  = Gpterm - Htr(i)/Ttr(i)*(T(iPT)-Ttr(i));
        end
        if tr_id == 1 || dPdTtr2 ~= 0
            if dPdTtr ~= 0
                Tc = Ttr1 + (P(iPT) - Pr)/dPdTtr;
            else
                Tc = Ttr1;
            end
            PtranT = Pr + (T(iPT) - Ttr1).*dPdTtr;
            if T(iPT)<Tc
                V(iPT)    = v0r;
                VdP(iPT)  = V(iPT)*(P(iPT) - Pr) + Vtr1*(PtranT - Pr);
                Gpterm    = 0;
            else
                V(iPT) = v0r + Vtr1;
            end
        end
        VdP(iPT) = VdP(iPT) + Gpterm;
    end
    if itype~=0%strcmp(phase(ip),'QUARTZ')==1 || strcmp(phase(ip),'COESITE')==1
        qphase    = ones(size(T))*2;
        Pstar     = zeros(size(T));
        Sstar     = zeros(size(T));
        V         = zeros(size(T));
        aa        = 0.549824e3; ba     = 0.65995;  ca     = -0.4973e-4; VPrTra = 22.688; Vdiff  = 2.047; k      = 38.5;  VPtTta = 0.23348e2; VPrTtb = 0.2372e2;    Stran  = 0.342;
        % set qphase = phase region of quartz
        qphase(T <= Ttr1 | P >= (Pr + k*(T-Ttr1))) = 1;
        % set Pstar and Sstar
        Pstar(T<=Ttr1) = Pr;
        Pstar(T>Ttr1&qphase==2) = P(T>Ttr1&qphase==2);
        Pstar(T>Ttr1&qphase==1) = Pr + k*(T(T>Ttr1&qphase==1)-Ttr1);
        V(qphase==2) = VPrTtb;
        V(qphase==1) = VPrTra + ca*(P(qphase==1)-Pr) + (VPtTta - VPrTra - ca*(P(qphase==1)-Pr)).*(T(qphase==1)-Tr) ./ ...
            (Ttr1 + (P(qphase==1)-Pr)/k - Tr);
        if itype == 2%strcmp(phase(ip),'COESITE')
            V = V - Vdiff;
        end
        % leading constant for [G,S]Vterm below is a conversion factor (cal/cm**3/bar)
        if itype == 1%strcmp(phase(ip),'QUARTZ')
            GVterm = 0.23901488e-1 * (VPrTra*(P-Pstar) + VPrTtb*(Pstar-Pr)  ...
                -         0.5*ca*(2*Pr*(P-Pstar) - ((P).^2-Pstar.^2)) ...
                -        ca*k*(T-Tr).*(P-Pstar) ...
                +        k*(ba + aa*ca*k)*(T-Tr).*log((aa + P/k)./(aa + Pstar/k)));
        else
            GVterm = 0.23901488e-1 * ((VPrTra-Vdiff)*(P-Pstar) ...
                +        (VPrTtb-Vdiff)*(Pstar-Pr) - 0.5*ca*(2*Pr*(P-Pstar) ...
                -        ((P).^2-Pstar.^2)) - ca*k*(T-Tr).*(P-Pstar) ...
                +        k*(ba + aa*ca*k)*(T-Tr).*log((aa + P/k)./(aa + Pstar/k)));
        end
        VdP = GVterm'*1/0.23901488;
    end
    VdP = VdP';
elseif eos == 11
    P   = P*1e3;Pr = 1;
    Vr  = td(8);
    VdP = Vr*(P-Pr);
    V   = Vr;
else
    V   = 0;
    VdP = 0;
end
end
function [int_VdP,V] = HP_CORK(Tarr,Parr,i_fl)
% CORK Equation of state and VdP from HP 1991 with Virial terms from HP 1998
R  = 8.314/1e3;
a_CORK    = [1113.4     659.8    ;
    -0.88517   0.21078   ;
    4.5300e-3 -6.3976e-4 ;
    -1.3183e-5  0;
    -0.22291    0;
    -3.8022e-4  0;
    1.7791e-7   0;
    5.8487      0;
    -2.1370e-2  0;
    6.8133e-5  0];
b_CORK     = [1.465 3.057];
vir_CORK   = [1.9853e-3; -8.9090e-2; 8.0331e-2];
vir_CORK_CO2 = [5.40776e-3  -1.59046e-6; % a virial
    -1.78198e-1 2.45317e-5   % b virial
    0          0]; % c virial
P0     = [2, 5];
V       = zeros(length(Tarr),1);
int_VdP = zeros(length(Tarr),1);
if i_fl == 1 || i_fl == 2 % for H2O or CO2
    for iT = 1:length(Tarr)
        T = Tarr(iT);
        P = Parr(iT);
        b = b_CORK(i_fl);
        if i_fl == 1
            dT = T-673;
            a_gas       = a_CORK(1,i_fl) + a_CORK(8:10,i_fl)'*(-dT).^[1 2 3]';
            Psat        = -13.627e-3 + 7.29395e-7*T^2 - 2.34622e-9*T^3 + 4.83607e-15*T^5;
            if P > Psat && T < 673
                a = a_CORK(1,i_fl) + a_CORK(2:4 ,i_fl)'*(-dT).^[1 2 3]';
            elseif P <= Psat && T < 673
                a = a_gas;
            else
                a = a_CORK(1,i_fl) + a_CORK(5:7 ,i_fl)'*( dT).^[1 2 3]';
            end
        elseif i_fl == 2 % CO2, see HP 1998 paper
            Psat        = 0;
            a           = a_CORK(1,i_fl) + a_CORK(2,i_fl)*T + a_CORK(3,i_fl)*T^2;
            vir_CORK    = vir_CORK_CO2*[1; T];
        end
        Vmrk = Vmrk_fun(T,P,a,b,R);
        Vvir = vir_CORK'*(P-P0(i_fl)).^[1;1/2;1/4];
        if T <= 695 && P > Psat && i_fl == 1
            Vmrk = min(Vmrk);
        else
            Vmrk = max(Vmrk);
        end
        if P > P0(i_fl)
            V(iT) = Vmrk + Vvir;
        else
            V(iT) = Vmrk;
        end
        %% For case A and B in appendix of HP 1991
        ln_g = ln_gam(T,P,Vmrk,a,b,vir_CORK(1),vir_CORK(2),vir_CORK(3),P0(i_fl),R); % has to be Vmrk only!
        %% For case C
        if i_fl == 1
            if T<695 && P>Psat
                if T<673
                    Vgas  = max(Vmrk_fun(T,Psat,a_gas,b,R));
                    ln_g1 = ln_gam(T,Psat,Vgas,a_gas,b,vir_CORK(1),vir_CORK(2),vir_CORK(3),P0(i_fl),R);
                else
                    Vgas  = max(Vmrk_fun(T,Psat,a,b,R));
                    ln_g1 = ln_gam(T,Psat,Vgas,a,b,vir_CORK(1),vir_CORK(2),vir_CORK(3),P0(i_fl),R);
                end
                Vliq  = min(Vmrk_fun(T,Psat,a,b,R));
                
                ln_g2 = ln_gam(T,Psat,Vliq,a,b,vir_CORK(1),vir_CORK(2),vir_CORK(3),P0(i_fl),R);
                ln_g  = ln_g1 - ln_g2 + ln_g;
            end
        end
        f_h2o       = exp(ln_g+log(P*1e3)); %%%% NOT SO OBVIOUS that this is P in bar
        int_VdP(iT) = R*T*log(f_h2o);
    end
elseif i_fl == 5 % O2
    int_VdP = 0;
elseif i_fl == 3 || i_fl == 4 || i_fl == 6
    T  = Tarr;
    P  = Parr;
    a0 = 5.45963e-5;
    a1 = -8.63920e-6;
    b0 = 9.18301e-4;
    c0 = -3.30558e-5;
    c1 = 2.30524e-6;
    d0 = 6.93054e-7;
    d1 = -8.38293e-8;
    if i_fl == 3 % CO (not implemented in THERMOCALC or Perple_X
        Tc = 132.9;
        Pc = 0.035;
    end
    if i_fl == 4 % CH4
        Tc = 190.6;
        Pc = 0.0460;
    end
    if i_fl == 6 % H2
        Tc = 41.2;
        Pc = 0.0211;
    end
    a = a0*Tc^(5/2)/Pc + a1*Tc^(3/2)/Pc*T;
    b = b0*Tc/Pc;
    c = c0*Tc/(Pc^(3/2)) + c1/Pc^(3/2)*T;
    d = d0*Tc/Pc^2 + d1/Pc^2*T;
    int_VdP = R*T.*log(1000*P) + b.*P + a./(b.*sqrt(T)).*( log(R*T + b.*P) - log(R*T + 2*b.*P) ) + 2/3*c.*P.*sqrt(P) + d/2.*P.^2;
end
end
function ln_g = ln_gam(T,P,V,a,b,a_vir,b_vir,c_vir,P0,R)
z = (P*V)/(R*T);
B = (b*P)/(R*T);
A = a/(b*R*T^(3/2));
ln_g_vir = (1/(R*T))*((4/5)*c_vir*(P-P0)^(5/4) + (2/3)*b_vir*(P-P0)^(3/2)+a_vir*((1/2)*(P-P0)^2));
ln_g     = z-1-log(z-B)-A*log(1+B/z);
if P>P0
    ln_g = ln_g + ln_g_vir;
end
end
function Vmrk = Vmrk_fun(T,P,a,b,R)
p_3 = P;
p_2 = -R*T;
p_1 = -(b*R*T+b^2*P-a/sqrt(T));
p_0 = -a*b/sqrt(T);

sol = roots([p_3 p_2 p_1 p_0]);
ind = zeros(size(sol));
for i = 1:length(sol)
    if isreal(sol(i))
        ind(i) = isreal(sol(i));
    end
end
Vmrk = sol(logical(ind));
end