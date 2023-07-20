clear,addpath ../
T     = 300 + 273.15;
P     = 1e9;
phase = {'Olivine','Orthopyroxene'};
Cname = {'Si','Mg','Fe','O'};
td    = init_thermo(phase,Cname,'solution_models_H18');
td(1).dz(:) = 1/2;
td(2).dz(:) = 1/2;
p     = props_generate(td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phase,td,p);