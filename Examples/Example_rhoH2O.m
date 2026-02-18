clear,addpath ../ ../Solutions ../EOS
T     = linspace(25,700,140) + 273.15;
P     = linspace(1e5,1e8,141);
phase = {'H2O,tc-ds633'};
[Temp,Pres] = ndgrid(T,P);
g     = tl_gibbs_energy(Temp(:),Pres(:),phase);
[rho_w,eps_di] = water_props(Temp(:),Pres(:),phase);
rho_IAPWS = water_props(Temp(:),Pres(:),phase,'IAPWS','S14');
td = init_thermo(phase);
g0     = tl_g0(Temp(:),Pres(:),td,rho_w,eps_di);
subplot(211),colormap jet,contourf(Temp-273.15,Pres/1e9,reshape(rho_w,length(T),length(P)),30);colorbar,shading interp
xlabel('T(\circC)'),ylabel('P(GPa)'),title('density of H_2O')
subplot(212),colormap jet,contourf(Temp-273.15,Pres/1e9,reshape(rho_IAPWS,length(T),length(P)),30);colorbar,shading interp
xlabel('T(\circC)'),ylabel('P(GPa)'),title('density of H_2O (IAPWS)')