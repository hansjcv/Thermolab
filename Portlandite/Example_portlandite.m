clear
T = linspace(298.15,1000,100);
P     = 1e5;
phs_name = {'Portlandite,Thermoddem'};
Cname = {'Ca','H','O'};
td    = init_thermo(phs_name,Cname);
[T2d,P2d] = ndgrid(T,P);
[g,Npc] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td);
S = -diff(g)./diff(T);
Tc = (T(1:end-1)+T(2:end))/2;
Cp = T(2:end-1).*diff(S)./diff(Tc);

a = 89.26;
b = 33.11e-3;
c = -10.36e5;
Cp_Thermoddem = a + b*T + c*T.^(-2);
plot(T(2:end-1),Cp,'d',T,Cp_Thermoddem,'o-')


% A = -284.929255;
% B = -0.044711;
% C = 21380.115459;
% D = 104.204551;
% E = -754249.169405;
logKd_Thermoddem = A + B*T + C*T.^(-1) + D*log10(T) + E*T.^(-2);