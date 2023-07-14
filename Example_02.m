clear, addpath ../
Temp   = linspace(300,1000,30) + 273.15;
Pres   = linspace(0.1,4,31)*1e9;
phase  = {'q,tc-ds633','coe,tc-ds633'};
[T,P]  = ndgrid(Temp,Pres);
g      = tl_gibbs_energy(T(:),P(:),phase);