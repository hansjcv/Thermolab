clear
T = linspace(298.15,550,100);
P     = 1e5;
R     = 8.3144;
phs_name = {'Portlandite,Thermoddem','Ca+2,supcrt','H+,supcrt','H2O,l,NIST'};
Cname = {'Ca','H','O','e'};
td    = init_thermo(phs_name,Cname);
[T2d,P2d] = ndgrid(T,P);
p             = props_generate(td);
[rho_w,eps_w] = water_props(T2d(:),P2d(:),phs_name,'JN91','JN91');
[g0,v0]       = tl_g0(T2d(:),P2d(:),td,rho_w,eps_w);
[g,Npc] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td,p,g0,v0,rho_w,eps_w);
A = -284.929255;
B = -0.044711;
C = 21380.115459;
D = 104.204551;
E = -754249.169405;
logKd_Thermoddem = A + B*T + C*T.^(-1) + D*log10(T) + E*T.^(-2);
v  = null(Npc);
v  = -v/v(1);
dg = v'*g;
Kd = exp(-dg./R./T);
plot(T,log10(Kd),T,logKd_Thermoddem,'o')