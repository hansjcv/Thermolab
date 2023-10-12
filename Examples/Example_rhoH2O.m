clear,addpath ../ ../Solutions ../EOS
T     = (25:10:1000) + 273.15;
P     = (0.001:0.001:0.1)*1e9;
phase = {'H2O,tc-ds633'};
[Temp,Pres] = ndgrid(T,P);
g     = tl_gibbs_energy(Temp(:),Pres(:),phase);
[rho_w,eps_di] = water_props(Temp(:),Pres(:),phase);
td = init_thermo(phase);
g0     = tl_g0(Temp(:),Pres(:),td,rho_w,eps_di);
colormap gray,contourf(Temp-273.15,Pres/1e9,reshape(rho_w,length(T),length(P)),30);colorbar,shading interp
xlabel('T(\circC)'),ylabel('P(GPa)'),title('density of H_2O')