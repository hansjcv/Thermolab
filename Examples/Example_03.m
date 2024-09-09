clear,addpath ../ ../Solutions/ ../EOS
Temp  = linspace(300,1000,30) + 273.15;
Pres  = linspace(0.1,4,31)*1e9;
phase = {'Olivine','Orthopyroxene','q,tc-ds633','per,tc-ds633'};
Cname = {'Si','Mg','Fe','O'};
td    = init_thermo(phase,Cname,'solution_models_H18');
[T,P] = ndgrid(Temp,Pres);
[g,Npc,pc_id] = tl_gibbs_energy(T(:),P(:),phase,td);