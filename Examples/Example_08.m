clear,clf, addpath ../ ../EOS ../Utilities/ ../Solutions/
T = 500 + 273.15;
P = 20e8;
phs_name      = {'Fluid-CO2-H2O'};
Cname         = {'Si','Al','Ca','Mg','Fe','Na','K','C','H','O'};
solname       = 'solution_models_HP98';
td            = init_thermo(phs_name,Cname,solname);
p             = props_generate(td);
[g0,v0]       = tl_g0(T,P,td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phs_name,td);
for ip = 1:length(phs_name)
    mu = tl_chemical_potential(T,P,td(ip),p{ip},g0(ip),v0(ip));
    chk(ip) = max(abs(g(pc_id==ip)-sum(p{ip}.*mu,2)));
end