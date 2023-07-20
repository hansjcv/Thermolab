clear,addpath ../ ../Solutions
T     = (300:100:1000) + 273.15;
P     = (0.1:0.25:1)*1e9;
phase = {'aqSi,tc-ds55','SiO2,aq,DEW'};
[Temp,Pres] = ndgrid(T,P);
g     = tl_gibbs_energy(Temp(:),Pres(:),phase);
[rho_w,eps_di] = water_props(Temp(:),Pres(:),phase);
td = init_thermo(phase);
g0     = tl_g0(Temp(:),Pres(:),td,rho_w,eps_di);