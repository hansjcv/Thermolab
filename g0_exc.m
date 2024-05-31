function dG = g0_exc(T,P,td,dgex,rho_w,eps_di)
if dgex == 1
    Tr = 25 + 273.15; Pr = 0.001;P = P/1e8;
    k298 = td(10);a0=td(8);
    Tc0  = td(14);Smax = td(15);Vmax   = td(16);
    %% Landau with mistake (see technical note)
    k_T       = k298*(1-1.5e-4*(T-Tr));
    Tc      = Tc0 + Vmax./Smax.*(P-1*Pr);
    Q_298   = (1-Tr./Tc0).^(1/4);
    Q       = (1-(T./Tc)).^(1/4);
    h_298   = Smax.*Tc0.*(Q_298.^2-(1/3)*Q_298.^6);
    s_298   = Smax.*Q_298.^2;
    v_t     = Vmax.*Q_298.^2.*(1+a0.*(T-Tr)-20*a0.*(sqrt(T)-sqrt(Tr)));
    intvdP  = (1/3)*v_t.*k_T.*((1+(4*(P-Pr))./k_T).^(3/4)-1);
    G_land  = Smax.*((T-Tc).*Q.^2 + (1/3).*Tc.*Q.^6);
    dG      = (h_298 - T.*s_298 + intvdP + G_land.*logical(T<=Tc))*1e3;
elseif dgex == 2
    Tr = 25 + 273.15; P = P/1e8;    
    %% Load data
    Tc0  = td(12);Smax = td(13);Vmax   = td(14);
    sfdh = td(15);sfdv = td(16);sfw    = td(17);sfwv = td(18);sfn = td(19);sffac = td(20); od = td(26);
    %% Landau corrected (see technical note)
    q2  = zeros(size(T,1),size(P,2));
    if ~isnan(Smax)
        q20       = sqrt(1-Tr./Tc0);
        tcp       = Tc0 + P.*(Vmax./Smax);
        temp_q2   = sqrt((tcp-T)./Tc0);
        q2(T<tcp) = temp_q2(T<tcp);
        g_exc     = Smax.*(Tc0.*(q20.*(1 - 1/3.*q20.^2) + 1/3*q2.^3) - q2.*tcp) - T.*Smax.*(q20 - q2) + P.*Vmax.*q20;
        if od == 1 % in case of ordered endmember ?
            g_exc     = Tc0.*Smax.*(-2/3 + q20.*(1 - q20.^2/3)) - Tc0.*Smax.*(q20 - 1) + P.*Vmax.*(q20 - 1);
        elseif od == 2 % in case of disordered endmember ?
            g_exc     = Tc0.*Smax.*q20.*(1 - q20.^2/3) - T.*Smax.*q20 + P.*Vmax.*q20;
        end
    else
        g_exc = 0;
    end
    %% Symmetric formalism (see Holland and Powell 1996)
    if ~isnan(sfn)
        R = 8.3144/1e3;
        if od == 0 % equilibrium
            for iT = 1:length(T(:))
                q(iT,1) = symm_form(T(iT),sfdh+P(iT)*sfdv,sfw + P(iT)*sfwv,sfn,sffac);
            end
            h_sf = sfdh + P.*sfdv + q.*(sfw-sfdh+P*(sfwv-sfdv))-q.^2.*(sfw+P*sfwv); % (1-Q)*sfdh + Q*(1-Q)*sfw
            s_sf = -sffac*R*((1+sfn*q).*log((1+sfn*q)/(sfn+1))+sfn*(1-q).*log((1-q)/(sfn+1)) ...
                +  sfn*(1-q).*log(sfn*(1-q)/(sfn+1))+sfn*(sfn+q).*log((sfn+q)/(sfn+1)) )/(sfn+1);
            g_sf = h_sf - T.*s_sf;
        elseif od == 2 % disordered
            if sffac<0
                s_sf = sffac * R * (log(1/(sfn+1)) + sfn*log(sfn/(sfn+1)))*(1/sffac-sfn)/(sfn+1);
            else
                s_sf = sffac * R * (log(1/(sfn+1)) + sfn*log(sfn/(sfn+1)));
            end
            g_sf = sfdh + T.*s_sf + P.*sfdv;
        elseif od == 1 % ordered
            g_sf = 0;
        end
    else
        g_sf = 0;
    end
    dG = (g_sf + g_exc)*1e3;
elseif dgex == 3
    Tref   = 298.15;P = P/1e8;
    eta    = 0.166027e6;ZPrTr = -0.1278034682d-1;YPrTr = -0.5798650444d-4;    
    td(10) = td(10)*1e5;
    a_g = [-2.037662, 5.747000e-3 -6.557892e-6];
    b_g = [6.107361 -1.074377e-2  1.268348e-5];
    f_ag = [0.3666666e2 -0.1504956e-9 0.5017997e-13];% f_ag1 = 3.666666e-16; in Table 4 from Shock et al 1992, but seems not to matter much
    gref = 0;
    P = P*1e3; TC = T - 273.15; rho_w_gcm3 = rho_w/1e3;
    Z = -1./eps_di;                                                                     %eq. 32 Johnson et al 1992
    % Computes g using equations given by Shock et al. (1992)
    a_g    = a_g(1) + a_g(2)*TC.^1 + a_g(3)*TC.^2;                                      %eq. 50 Johnson et al 1992
    b_g    = b_g(1) + b_g(2)*TC.^1 + b_g(3)*TC.^2;                                      %eq. 51 Johnson et al 1992
    f      = (((TC - 155)/300).^(4.8) + f_ag(1)*((TC - 155)/300).^(16 )).* ...
        (f_ag(2)*(1000 - P).^3 + f_ag(3)*(1000 - P).^4); f((TC < 155) | (P > 1000) |  (TC > 355)) = 0;      %eq. 52 Johnson et al 1992
    g_shok = a_g.*(1 - rho_w_gcm3).^b_g; g_shok(rho_w_gcm3 >= 1) = 0;                   %eq. 49 Johnson et al 1992
    g_shok = g_shok  - f;                                                               %eq. 49 Johnson et al 1992
    % Computes the conventional Born coefficient (w) of the current aqueous species as a function of g, wref, and charge.
    if td(11)==0 || sum(td(1:10))==0 % neutral aqueous species or H+
        w     = td(10);
    else %        charged aqueous species other than H+
        reref = td(11)^2 / (td(10)/eta + td(11)/(3.082 + gref));%eq. 56 Johnson et al 1992
        re    = reref + abs(td(11)) * g_shok;                                   %eq. 48 Johnson et al 1992
        w     = eta * (td(11)^2./re - td(11)./(3.082 + g_shok));        %eq. 55 Johnson et al 1992
    end
    % HKF standard molal gibbs energies
    dG =  (w.*(-Z - 1) ...
        - td(10)*(-ZPrTr - 1) ...
        + td(10)*YPrTr*(T-Tref))*4.184;                                                %from eq. 59 Johnson et al 1992 see Dolejs RimG 76
elseif dgex == 4 || dgex == 5
    P = P/1e8;
    rho_gcm3 = rho_w*1e-3;
    Vref = 1.817884158651069;    
    Mw   = 18.0150;
    rho0 = Mw/10/Vref;
    Href = td(1);Sref  = td(2 );V1_298 = td(3); b = td(5 );
    if dgex == 4
        Cp0 = td(14);
    elseif dgex == 5
        Cp0 = td(21);
    end
    Tref = 298.15; Pref = 1e-3; T1 = T; T1(T>500) = 500; dadT = 9.5714e-6; a0 = 25.93e-5; B0 = 45.23e-6*1e3; % given in bar-1 convert to kbar-1    
    % Gibbs calculation aqueous species HP 1998
    Cps   = Cp0 - Tref*b;
    Gref  = Href - Tref*Sref;
    dG  = ((Href - T.*Sref + (P-Pref).*V1_298 + b.*(Tref.*T - 0.5*Tref^2 - 0.5*T.^2)...
        +  (Cps./(Tref.*dadT)).*(a0.*(T-Tref) - B0.*(P-Pref) + (T./T1).*log(rho_gcm3./rho0)) - Gref))*1e3;
elseif dgex == 6
    P = P/1e5;
    T1_lamb = td(16);Tref_lamb = td(17);l1 = td(18);l2 = td(19);dTdP = td(20);
    TP_lamb = T1_lamb + dTdP*(P-1);
    T_d = T1_lamb-TP_lamb;
    tr = Tref_lamb-T_d;    
    x1 = l1^2*T_d + 2*l1*l2*T_d.^2 + l2^2*T_d.^3;
    x2 = l1^2 + 4*l1*l2*T_d + 3*l2^2*T_d.^2;
    x3 = 2*l1*l2 + 3*l2^2*T_d;
    x4 = l2^2;
    Hlamb = x1.*(T-tr) + x2./2.*(T.^2-tr.^2) + x3./3.*(T.^3-tr.^3) + x4./4.*(T.^4-tr.^4);
    Slamb = x1.*(log(T)-log(tr)) + x2.*(T-tr) + x3./2.*(T.^2-tr.^2) + x4./3.*(T.^3-tr.^3);
    Glamb = Hlamb - T.*Slamb;
    Hlamb2 = x1.*(TP_lamb-tr)+x2./2.*(TP_lamb.^2-tr.^2) + x3./3.*(TP_lamb.^3-tr.^3) + x4./4.*(TP_lamb.^4-tr.^4);
    Slamb2 = x1.*(log(TP_lamb)-log(tr)) + x2.*(TP_lamb-tr) + x3./2.*(TP_lamb.^2-tr.^2) + x4./3.*(TP_lamb.^3-tr.^3);
    Glamb2 = Hlamb2 - T.*Slamb2;
    if T1_lamb~=0
        dG =  Glamb.*logical(T<=TP_lamb) + Glamb2.*logical(T>TP_lamb);
        if dTdP==0
            dG = - (T-TP_lamb).*Slamb.*logical(T<=TP_lamb);
        end
    end
elseif dgex == 7
    P = P/1e5; Pr = 1;
    d0 = td(21);d1 = td(22);d2 = td(23);d3 = td(24);d4 = td(25);d5 = td(26);Tmin = td(27);Tmax = td(28);
    Hds   = d0*(T-Tmin) + 2*d1*(T.^(0.5)-Tmin^0.5) - d2*(T.^(-1)-Tmin^(-1)) + d3*(log(T)-log(Tmin)) + d4*(T.^2-Tmin^2)/2 + d5*(T.^3-Tmin^3)/3;
    Sds   = d0*(log(T)-log(Tmin)) - 2*d1*(T.^(-0.5)-Tmin^(-0.5)) - d2*(T.^(-2)-Tmin^(-2))/2 + d3*(T-Tmin) + d4*(T-Tmin) + d5*(T.^2-Tmin^2)/2;
    Vds   = 0;
    dGds  = Hds - T.*Sds + Vds.*(P-Pr);
    Hds2  = d0*(T-Tmax) + 2*d1*(T.^(0.5)-Tmax^0.5) - d2*(T.^(-1)-Tmax^(-1)) + d3*(log(T)-log(Tmax)) + d4*(T.^2-Tmax^2)/2 + d5*(T.^3-Tmax^3)/3;
    Sds2  = d0*(log(T)-log(Tmax)) - 2*d1*(T.^(-0.5)-Tmax^(-0.5)) - d2*(T.^(-2)-Tmax^(-2))/2 + d3*(T-Tmax) + d4*(T-Tmax) + d5*(T.^2-Tmax^2)/2;   
    Vds2  = 0;   
    dGds2 = Hds2 - T.*Sds2 + Vds2.*(P-Pr);
    dG    =  dGds  + dGds2.*logical(T>Tmax);    
else
    dG = zeros(size(T,1),size(P,2));
end
end
function f = bragg_will(dH,W,Q,fac,n,T)
R = 8.3144/1e3;
% From paper:
K = ((n-n*Q).*(1-Q))./((1+n*Q).*(n+Q));
f = dH + W.*(2.*Q-1) + fac.*n./(n+1).*R.*T.*log(K);
end
function ihq = symm_form(T,dH,W,n,fac)

small = eps;
itlim = 50;

q1    = small;
q2    = 1-small;

f1 = bragg_will(dH,W,q1,fac,n,T);
f2 = bragg_will(dH,W,q2,fac,n,T);

if f2>0
    ihq = 1-small;
elseif f1<0
    ihq = small;
else
    qk = 0.5;
    fk = bragg_will(dH,W,qk,fac,n,T);
    for it = 1:itlim
        if fk>0
            q1 = qk;
        else
            q2 = qk;
        end
        qk = (q1+q2)/2;
        fk = bragg_will(dH,W,qk,fac,n,T);
        ihq = qk;
    end
end
end