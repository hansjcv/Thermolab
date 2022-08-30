function rho_w = rho_H2O(T,P,EOS)
P     = P/1e8;
rho_w = zeros(length(T),1);
switch EOS
    case 'JN91'
        %% Johnson and Norton 1991
        rho_g1 = 500;
        rho_g2 = 1000;
        rho_g3 = 1600;
%         for i = 1:length(T)
%             for j = 1:length(P)
%                 if P(j)<2
%                     if T(i) > 573.15
%                         rho_w(i,j) = fzero(@(rho) EOS_JN91(rho,T(i))-P(j),rho_g1); % T in K , P in kb
%                     else
%                         rho_w(i,j) = fzero(@(rho) EOS_JN91(rho,T(i))-P(j),rho_g2); % T in K , P in kb
%                     end
%                 else
%                     rho_w(i,j) = fzero(@(rho) EOS_JN91(rho,T(i))-P(j),rho_g3); % T in K , P in kb
%                 end
%             end            
%         end
        for i = 1:length(T)            
            if P(i)<2
                if T(i) > 573.15
                    rho_w(i) = fzero(@(rho) EOS_JN91(rho,T(i))-P(i),rho_g1); % T in K , P in kb
                else
                    rho_w(i) = fzero(@(rho) EOS_JN91(rho,T(i))-P(i),rho_g2); % T in K , P in kb
                end
            else
                rho_w(i) = fzero(@(rho) EOS_JN91(rho,T(i))-P(i),rho_g3); % T in K , P in kb
            end
        end
    case 'IAPWS'
        %% IAPWS EOS
        rho_g1 = 800;
        rho_g2 = 1100;
        rho_g3 = 1600;
%         for i = 1:length(T)
%             for j = 1:length(P)
%                 if P(j)<2
%                     if T(i) > 573.15
%                         rho_w(i,j) = fzero(@(rho) rho_IAPWS(rho,T(i),P(j)*1e2),rho_g1); % T in K , P in MPa
%                     else
%                         rho_w(i,j) = fzero(@(rho) rho_IAPWS(rho,T(i),P(j)*1e2),rho_g2); % T in K , P in MPa
%                     end
%                 else
%                     rho_w(i,j) = fzero(@(rho) rho_IAPWS(rho,T(i),P(j)*1e2),rho_g3); % T in K , P in MPa
%                 end
%             end
%         end
        for i = 1:length(T)           
            if P(i)<2
                if T(i) > 573.15
                    rho_w(i) = fzero(@(rho) rho_IAPWS(rho,T(i),P(i)*1e2),rho_g1); % T in K , P in MPa
                else
                    rho_w(i) = fzero(@(rho) rho_IAPWS(rho,T(i),P(i)*1e2),rho_g2); % T in K , P in MPa
                end
            else
                rho_w(i) = fzero(@(rho) rho_IAPWS(rho,T(i),P(i)*1e2),rho_g3); % T in K , P in MPa
            end
        end
    case 'ZD09'
        %% Zhang and Duan 2009 EOS
        rho_g1 = 1800;
        rho_g2 = 1000;
%         for i = 1:length(T)
%             for j = 1:length(P)
%                 if T(i) > 573.15
%                     rho_w(i,j) = fzero(@(rho) rho_ZD09(rho,T(i),P(j)),rho_g1); % T in K , P in kbar
%                 else
%                     warning('Out of calibration range for Zhang and Duan EOS, switching to IAPWS');
%                     if P(j)<2
%                         rho_w(i,j) = fzero(@(rho) rho_IAPWS(rho,T(i),P(j)*1e2),rho_g2); % T in K , P in MPa
%                     else
%                         rho_w(i,j) = fzero(@(rho) rho_IAPWS(rho,T(i),P(j)*1e2),rho_g1); % T in K , P in MPa
%                     end
%                 end
%             end
%         end
        for i = 1:length(T)
            if T(i) > 573.15
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
%         for i = 1:length(T)
%             for j = 1:length(P)
%                 if T(i) > 573.15
%                     rho_w(i,j) = fzero(@(rho) rho_ZD05(rho,T(i),P(j)),rho_g1); % T in K , P in kbar
%                 else
%                     rho_w(i,j) = fzero(@(rho) rho_ZD05(rho,T(i),P(j)),rho_g2); % T in K , P in kbar               
%                 end                
%             end
%         end
        for i = 1:length(T)            
            if T(i) > 573.15
                rho_w(i) = fzero(@(rho) rho_ZD05(rho,T(i),P(i)),rho_g1); % T in K , P in kbar
            else
                rho_w(i) = fzero(@(rho) rho_ZD05(rho,T(i),P(i)),rho_g2); % T in K , P in kbar
            end
        end
    case 'CORK'
        %% Holland and Powell (1991) CORK
        td = zeros(1,18);
        td(17) = 1;
%         [T2d,P2d] = ndgrid(T,P);
%         [~,V] = intVdP(T2d(:),P2d(:)*1e8,td,4);
%         rho_w = (18.0150/10)./V;       
%         rho_w = reshape(rho_w,length(T),length(P));        
        [~,V] = intVdP(T,P*1e8,td,4);
        rho_w = (18.0150e3/10)./V;               
end
end
function Pr_JN91  = EOS_JN91(rho,T)
% Data for Johnson and Norton 1991 H2O Eos
% Constants
R = 0.46152; % bar cm^3/g/K
gam1 = 13;
n    = 36;
kr   = 1;
T0   = 647.074;% K
% Fitting parameters Pbase
b  = [0.7478629,-0.3540782, 0        ,0.007159876, 0           ,-0.003528426];%cm3/g
B  = [1.1278334,-0.5944001, -5.010996,0          , 0.63684256  ,0           ];%cm3/g
% Fitting parameters Presid
kk = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 9 9 9 9 3 3 1 5 2 2 2 4];
ll = [1 2 4 6 1 2 4 6 1 2 4 6 1 2 4 6 1 2 4 6 1 2 4 6 1 2 4 6 1 2 4 6 0 3 3 3 0 2 0 0];                     
g  = [-.53062968529023d3,  .22744901424408d4,  .78779333020687d3 ...
    -.69830527374994d2,  .17863832875422d5, -.39514731563338d5 ...
    .33803884280753d5, -.13855050202703d5, -.25637436613260d6 ...
    .48212575981415d6, -.34183016969660d6,  .12223156417448d6 ...
    .11797433655832d7, -.21734810110373d7,  .10829952168620d7 ...
    -.25441998064049d6, -.31377774947767d7,  .52911910757704d7 ...
    -.13802577177877d7, -.25109914369001d6,  .46561826115608d7 ...
    -.72752773275387d7,  .41774246148294d6,  .14016358244614d7 ...
    -.31555231392127d7,  .47929666384584d7,  .40912664781209d6 ...
    -.13626369388386d7,  .69625220862664d6, -.10834900096447d7 ...
    -.22722827401688d6,  .38365486000660d6,  .68833257944332d4 ...
    .21757245522644d5, -.26627944829770d4, -.70730418082074d5 ...
    -.22500000000000d0, -.16800000000000d1 ...
    .5500000000000d-1, -.93000000000000d2];
% Fitting parameters Presid in critical region
T_i   = [zeros(1,n) .64d3,  .64d3,  .6416d3, .27d3];
rho_i = [zeros(1,n) .319d0, .319d0, .319d0,  .155d1];
beta  = [zeros(1,n) .2d5,   .2d5,   .4d5,    .25d2];
alp_i = [zeros(1,n) .34d2,  .4d2,   .3d2,    .105d4];

rho  = rho*1e-3; 
b_v2 = b(2)*log(T/T0);
B_v2 = 0;
for i = 1:length(b)
    j = i-1;
    if j~=1
        b_v2 = b_v2  + b(i)*(T0./T).^(j);
    end
    B_v2 = B_v2 + B(i)*(T0./T).^(j);
end
y     = b_v2.*rho/4;
Pbase = R*T.*rho.*((1 + (gam1 - 2).*y + (1 - gam1 + gam1.^2/3).*y.^2)./((1 - y).^3) ...
    +                4*y.*(B_v2./b_v2 - (1 + gam1)/4));
Presid = zeros(size(rho));
for i = 1:n
    Presid = Presid + g(i)*(T0./T).^ll(i).*rho.^2.*(1-exp(-kr.*rho)).^(kk(i)-1).*exp(-kr*rho);
end
for i = 37:40
    del = (rho-rho_i(i))/rho_i(i);
    tau = (T-T_i(i))/T_i(i);
    Presid = Presid + ...
        g(i)/rho_i(i)*rho.^2.*exp(-alp_i(i)*del.^kk(i)-beta(i)*tau.^2).* ...
        (ll(i)*del.^(ll(i)-1)-alp_i(i)*kk(i)*del.^(kk(i)+ll(i)-1));
end
Pr_JN91 = (Pbase + Presid)/100; % convert to kbar
end
function Pr_ZD05  = rho_ZD05(density,TK,Pin)
density = density/1e3; % convert to g/cm3;

% Zhang and Duan 2005 EOS for water
m = 18.01528;
ZD05_R = 83.14467;            %'Gas Constant in units of cm^3 bar/mol/K
ZD05_Vc = 55.9480373;        %'Critical volume in units of cm^3/mol
ZD05_Tc = 647.25;            %'Critical temperature in units of Kelvin

Vr = m ./ density / ZD05_Vc;
Tr = TK / ZD05_Tc;

B = 0.349824207 - 2.91046273 ./ (Tr .* Tr) + 2.00914688 ./ (Tr .* Tr .* Tr);
C = 0.112819964 + 0.748997714 ./ (Tr .* Tr) - 0.87320704 ./ (Tr .* Tr .* Tr);
D = 0.0170609505 - 0.0146355822 ./ (Tr .* Tr) + 0.0579768283 ./ (Tr .* Tr .* Tr);
E = -0.000841246372 + 0.00495186474 / (Tr .* Tr) - 0.00916248538 / (Tr .* Tr .* Tr);
f = -0.100358152 ./ Tr;
g = -0.00182674744 .* Tr;

delta = 1 + B ./ Vr + C ./ (Vr .* Vr) + D ./ Vr.^4 + E ./ Vr.^5 + (f ./ (Vr .* Vr) + g ./ Vr.^4) * exp(-0.0105999998 ./ (Vr .* Vr));


Pr_ZD05 = ZD05_R .* TK .* density .* delta / m / 1e3 - Pin;% converted to kbar
end
function Pr_ZD09  = rho_ZD09(rho,T,Pin)
temperature = T - 273.15;
m = 18.01528;
ZD09_R = 0.083145;        %'Gas constant in units of dm^3 bar/mol/K
ZD09_epsilon = 510;       %'Lenard-Jones parameter in units of K
ZD09_omega = 2.88;        %'Lenard-Jones parameter in units of 1E-10 m
ZD09_c1 = 3.0636*ZD09_omega^3/ZD09_epsilon;
density = rho*1e3/1e6; % g/cm3
%'Prefactor calculated from 1000 * pow(ZD09_omega / 3.691, 3)
dm = 475.05656886 * density;
%'Prefactor calculated from 0.001 * pow(3.691 / ZD09_omega, 3)
Vm = 0.0021050125 * (m ./ density);
%'Prefactor calculated from 154 / ZD09_epsilon
Tm = 0.3019607843 * (temperature + 273.15);  % 'temperature must be converted to Kelvin

B = 0.029517729893 - 6337.56452413 ./ (Tm .* Tm) - 275265.428882 ./ (Tm .* Tm .* Tm);
C = 0.00129128089283 - 145.797416153 ./ (Tm .* Tm) + 76593.8947237 ./ (Tm .* Tm .* Tm);
D = 2.58661493537E-06 + 0.52126532146 ./ (Tm .* Tm) - 139.839523753 ./ (Tm .* Tm .* Tm);
E = -2.36335007175E-08 + 0.00535026383543 ./ (Tm .* Tm) - 0.27110649951 ./ (Tm .* Tm .* Tm);
f = 25038.7836486 ./ (Tm .* Tm .* Tm);

delta = 1 + B ./ Vm + C ./ (Vm .* Vm) + D ./ Vm.^4 + E ./ Vm.^5 + ...
    f ./ (Vm .* Vm) .* (0.73226726041 + 0.015483335997 ./ (Vm .* Vm)) .* exp(-0.015483335997 ./ (Vm .* Vm));

Pm = ZD09_R * Tm .* delta ./ Vm;

Pr_ZD09 = Pm / ZD09_c1 / 1e3 - Pin; % kbar!
end
function [Pr_IAPWS,F] = rho_IAPWS(rho,T,Pin)
% Parameters
rho_c = 322; % kg/m3
Tc    = 647.096; % K
R     = 0.46151805; % kJ/kg/K
n0    = [-8.3204464837497  6.6832105275932 3.00632 0.012436   0.97315    1.27950    0.96956    0.24873   ];
gam0  = [nan               nan             nan     1.28728967 3.53734222 7.74073708 9.24437796 27.5075105];
eq6_pars = [
    1             nan            1             -0.5              0.12533547935523e-1 0 0 0 0 0 0 0 0 0 0
    2             nan            1            0.875              0.78957634722828e1 0 0 0 0 0 0 0 0 0 0
    3             nan            1               1              -0.87803203303561e1 0 0 0 0 0 0 0 0 0 0
    4             nan            2              0.5              0.31802509345418 0 0 0 0 0 0 0 0 0 0
    5             nan            2             0.75             -0.26145533859358 0 0 0 0 0 0 0 0 0 0
    6             nan            3            0.375             -0.78199751687981e-2 0 0 0 0 0 0 0 0 0 0
    7             nan            4               1               0.88089493102134e-2 0 0 0 0 0 0 0 0 0 0
    8             1              1               4              -0.66856572307965 0 0 0 0 0 0 0 0 0 0
    9             1              1               6               0.20433810950965 0 0 0 0 0 0 0 0 0 0
    10            1              1              12              -0.66212605039687e-4 0 0 0 0 0 0 0 0 0 0
    11            1              2               1              -0.19232721156002 0 0 0 0 0 0 0 0 0 0
    12            1              2               5              -0.25709043003438 0 0 0 0 0 0 0 0 0 0
    13            1              3               4               0.16074868486251 0 0 0 0 0 0 0 0 0 0
    14            1              4               2              -0.40092828925807e-1 0 0 0 0 0 0 0 0 0 0
    15            1              4              13               0.39343422603254e-6 0 0 0 0 0 0 0 0 0 0
    16            1              5               9              -0.75941377088144e-5 0 0 0 0 0 0 0 0 0 0
    17            1              7               3               0.56250979351888e-3 0 0 0 0 0 0 0 0 0 0
    18            1              9               4              -0.15608652257135e-4 0 0 0 0 0 0 0 0 0 0
    19            1             10             11                0.11537996422951e-8 0 0 0 0 0 0 0 0 0 0
    20            1             11              4                0.36582165144204e-6 0 0 0 0 0 0 0 0 0 0
    21            1             13             13               -0.13251180074668e-11 0 0 0 0 0 0 0 0 0 0
    22            1             15              1               -0.62639586912454e-9 0 0 0 0 0 0 0 0 0 0
    23            2              1               7              -0.10793600908932 0 0 0 0 0 0 0 0 0 0
    24            2              2               1               0.17611491008752e-1 0 0 0 0 0 0 0 0 0 0
    25            2              2               9               0.22132295167546 0 0 0 0 0 0 0 0 0 0
    26            2              2              10              -0.40247669763528 0 0 0 0 0 0 0 0 0 0
    27            2              3              10               0.58083399985759 0 0 0 0 0 0 0 0 0 0
    28            2              4               3               0.49969146990806e-2 0 0 0 0 0 0 0 0 0 0
    29            2              4               7              -0.31358700712549e-1 0 0 0 0 0 0 0 0 0 0
    30            2              4              10              -0.74315929710341 0 0 0 0 0 0 0 0 0 0
    31            2              5              10               0.47807329915480 0 0 0 0 0 0 0 0 0 0
    32            2              6               6               0.20527940895948e-1 0 0 0 0 0 0 0 0 0 0
    33            2              6              10              -0.13636435110343 0 0 0 0 0 0 0 0 0 0
    34            2              7              10               0.14180634400617e-1 0 0 0 0 0 0 0 0 0 0
    35            2              9               1               0.83326504880713e-2 0 0 0 0 0 0 0 0 0 0
    36            2              9               2              -0.29052336009585e-1 0 0 0 0 0 0 0 0 0 0
    37            2              9               3               0.38615085574206e-1 0 0 0 0 0 0 0 0 0 0
    38            2              9               4              -0.20393486513704e-1 0 0 0 0 0 0 0 0 0 0
    39            2              9               8              -0.16554050063734e-2 0 0 0 0 0 0 0 0 0 0
    40            2             10              6                0.19955571979541e-2 0 0 0 0 0 0 0 0 0 0
    41            2             10              9                0.15870308324157e-3 0 0 0 0 0 0 0 0 0 0
    42            2             12              8               -0.16388568342530e-4 0 0 0 0 0 0 0 0 0 0
    43            3              3              16               0.43613615723811e-1 0 0 0 0 0 0 0 0 0 0
    44            3              4              22               0.34994005463765e-1 0 0 0 0 0 0 0 0 0 0
    45            3              4              23              -0.76788197844621e-1 0 0 0 0 0 0 0 0 0 0
    46            3              5              23               0.22446277332006e-1 0 0 0 0 0 0 0 0 0 0
    47            4             14             10               -0.62689710414685e-4 0 0 0 0 0 0 0 0 0 0
    48            6              3              50              -0.55711118565645e-9 0 0 0 0 0 0 0 0 0 0
    49            6              6              44              -0.19905718354408 0 0 0 0 0 0 0 0 0 0
    50            6              6              46               0.31777497330738 0 0 0 0 0 0 0 0 0 0
    51            6              6               50             -0.11841182425981 0 0 0 0 0 0 0 0 0 0
    52            nan            3               0              -0.31306260323435e2 20 150 1.21 1 0 0 0 0 0 0
    53            nan            3               1              -0.31546140237781e2 20 150 1.21 1 0 0 0 0 0 0
    54            nan            3               4              -0.25213154341695e4 20 250 1.25 1 0 0 0 0 0 0
    55            0              0               0              -0.14874640856724   0  0.3   0    0  3.5 0.85 0.2 28 700 0.32
    56            0              0               0               0.31806110878444   0  0.3   0    0  3.5 0.95 0.2 32 800 0.32
    ];
c = eq6_pars(:,2);
d = eq6_pars(:,3);
t = eq6_pars(:,4);
n = eq6_pars(:,5);
alph  = eq6_pars(:,6);
beta = eq6_pars(:,7);
gam  = eq6_pars(:,8);
epsi = eq6_pars(:,9);
a  = eq6_pars(:,10);
b  = eq6_pars(:,11);
B  = eq6_pars(:,12);
C  = eq6_pars(:,13);
D  = eq6_pars(:,14);
A  = eq6_pars(:,15);
% Preprocess
delta = rho/rho_c;
tau   = Tc./T;
% Initialize
phi0  = 0;
phir = 0;
dphir = 0;
% The ideal part of Helmholtz energy
for i = 4:8
    phi0  = phi0 + n0(i)*log(1-exp(-gam0(i)*tau));
end
phi0 = phi0 + log(delta) + n0(1) + n0(2)*tau + n0(3)*log(tau);
% The residual part
for i = 1:7
    phir  = phir + n(i)*delta.^d(i).*tau.^t(i);
    dphir  = dphir + n(i)*d(i)*delta.^(d(i)-1).*tau.^t(i);
end
for i = 8:51
    phir  = phir  + n(i)*delta.^d(i).*tau.^t(i).*exp(-delta.^c(i));
    dphir = dphir + n(i)*exp(-delta.^c(i)).*(delta.^(d(i)-1).*tau.^t(i).*(d(i)-c(i).*delta.^c(i)));
end
for i = 52:54
    phir   = phir  + n(i)*delta.^d(i).*tau.^t(i).*exp(-alph(i).*(delta-epsi(i)).^2-beta(i)*(tau-gam(i)).^2);
    dphir  = dphir + n(i)*delta.^d(i).*tau.^t(i).*exp(-alph(i).*(delta-epsi(i)).^2-beta(i)*(tau-gam(i)).^2) ...
        .*(d(i)./delta-2*alph(i)*(delta-epsi(i)));
end
for i = 55:56
    theta = (1-tau) + A(i)*((delta-1).^2).^(1/(2*beta(i)));
    delta_L = theta.^2+B(i)*((delta-1).^2).^a(i);
    psi = exp(-C(i)*(delta-1).^2-D(i)*(tau-1).^2);
    phir   = phir + n(i)*delta_L.^b(i).*delta.*psi;
    dpsi_d = -2*C(i)*(delta-1).*psi;
    d_D_dd = (delta-1).*(A(i).*theta*2/beta(i).*((delta-1).^2).^(1/(2*beta(i))-1) + 2*B(i)*a(i)*((delta-1).^2).^(a(i)-1) );
    d_Dbi_d = b(i)*delta_L.^(b(i)-1).*d_D_dd;
    dphir  = dphir + n(i)*( delta_L.^b(i).*(psi + delta.*dpsi_d) + d_Dbi_d.*delta.*psi);
end
F = (phi0 + phir)*R*1e3.*T;
% Pressure
Pr_IAPWS    = (1 + delta.*dphir).*rho*R*1e3.*T/1e6 - Pin;
end