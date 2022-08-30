function [SdT,S] = intSdT(T,P,td,ceos)
if ceos == 1
    Sr = td(2);a = td(4);b = td(5);c = td(6);d = td(7);
    P  = P/1e8;
    Tr = 25 + 273.15;
    S   = Sr +              log(T/Tr)*a...
        +                      (T-Tr)*b ...
        + (1./(2*Tr^2) - 1./(2*T.^2))*c ...
        +  (2./sqrt(Tr) - 2./sqrt(T))*d;
    SdT =                      (T.*log(T/Tr) + Tr - T)*a ...
        +                              1/2*((T-Tr).^2)*b ...
        +                    ((T - Tr).^2./(2*T*Tr^2))*c ...
        +    (-4*sqrt(T) + 2*sqrt(Tr) + 2*T./sqrt(Tr))*d ...
        +                                       (T-Tr)*Sr;
    S   = S*1e3; % convert to Joules
    SdT = SdT*1e3;% convert to Joules
elseif ceos == 2
    Tref   = 298.15;
    theta = 228;
    td(9)  = td(9)*1e4;
    SdT =  (td(3)*(T-Tref) ...
        +  td(8)*(T.*log(T/Tref)-T+Tref) ...
        +  td(9)* ((1./(T-theta) - 1/(Tref-theta)).* ((theta-T)/theta) ...
        -                   T/theta^2.*log((Tref*(T-theta))./(T*(Tref-theta)))))*4.184;
    S =   (td(3)+ ...
        td(8)*log(T/Tref) ...
        - td(9)*(-(T - Tref)./((T - theta)*(Tref - theta)*theta) ...
        + log(Tref*(T - theta)./(T*(Tref - theta)))/theta^2))*4.184;
elseif ceos == 3 % SUPCRT Minerals
    Tr = 298.15;
    g0r = td(1);h0r = td(2);Sr = td(3); v0r = td(4);a = td(5);b = td(6);c = td(7);
    Ttr1     = td(8); Htr1  = td(9);Vtr1 = td(10); dPdTtr = td(11); atr = td(12); btr = td(13); ctr = td(14);
    Ttr2    = td(15);Htr2 = td(16);Vtr2 = td(17);dPdTtr2 = td(18);atr2 = td(19);btr2 = td(20);ctr2 = td(21);
    Ttr3    = td(22);Htr3 = td(23);Vtr3 = td(24);dPdTtr3 = td(25);atr3 = td(26);btr3 = td(27);ctr3    = td(28);
    Ttr4    = td(29);tr_id   = td(30);  itype   = td(31);
    Ttr = [Tr Ttr1 Ttr2 Ttr3 ];Htr = [0  Htr1 Htr2 Htr3 ];Vtr = [v0r Vtr1 Vtr2 Vtr3 ];atr = [a   atr atr2 atr3 ];btr = [b   btr btr2 btr3 ];ctr = [c   ctr ctr2 ctr3 ];
    for iT = 1:length(T)
        S(iT) = Sr;
        H_atP = 0;
        Tlo = Ttr(abs(atr)>0);
        Thi = [Tlo(2:end) T(iT)];
        for i = 1:tr_id+1
            T1 = Tlo(i); T2 = min(T(iT),Thi(i));
            if T1>T(iT), break,end
            S(iT)   = S(iT) + atr(i)*log(T2/T1) + btr(i)*(T2-T1)       - ctr(i)/2*(1/T2^2 - 1/T1^2);
            H_atP   = H_atP + atr(i)*(T2-T1)    + btr(i)/2*(T2^2-T1^2) - ctr(i)*(1/T2-1/T1);
        end
        SdT(iT,1) = - Sr*Tr - H_atP + T(iT)*S(iT);
    end
elseif ceos == 4 % Shomate
    t   = T/1000;
    H0  = td(2)*t + td(3)*t.^2/2 + td(4)*t.^3/3 + td(5)*t.^4/4 - td(6)./t + td(7) ;
    S   = td(2)*log(t) + td(3)*t + td(4)*t.^2/2 + td(5)*t.^3/3 - td(6)./(2*t.^2) + td(8) + 0.2;
    SdT = -(H0 - T.*S*1e-3 + 69.5584)*1e3;
elseif ceos == 5
    Tr = 298.15;
    Sr = td(2);a = td(3);b = td(4); c = td(5);d = td(6); e = td(7);
     SdT = (-T*log(Tr) + T.*log(T) + Tr - T)*a ...
         + (1/2*Tr^2 - Tr*T + 1/2*T.^2)*b ...
         + (1./(2*T) - 1/Tr + T/(2*Tr^2))*c ...
         + (1/3*Tr^3 - 1/2*Tr^2*T + 1/6*T.^3)*d ...
         + (-4*sqrt(T) + 2*sqrt(Tr) + 2*T/sqrt(Tr))*e - Tr*Sr + T*Sr;
     S = (log(T) - log(Tr))*a + (T - Tr)*b + (-1./(2*T.^2) + 1/(2*Tr^2))*c + (-Tr^2/2 + T.^2/2)*d + (-2./sqrt(T) + 2/sqrt(Tr))*e + Sr;
else
    S   = 0;
    SdT = 0;
end
end