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