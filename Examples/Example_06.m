clear,addpath ../ ../Solutions/
T     = 800 + 273.15;
P     = 1e9;
phase = {'Olivine','Orthopyroxene'};
Cname = {'Si','Mg','Fe','O'};
td    = init_thermo(phase,Cname,'solution_models_H18');
td(1).dz(:) = 1/2;
td(2).dz(:) = 1/2;
p             = props_generate(td);
[rho_w,eps_w] = water_props(T,P,phase,'PS94','S14')
[g0,v0]       = tl_g0(T,P,td,rho_w,eps_w);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phase,td,p,g0,v0,rho_w,eps_w);