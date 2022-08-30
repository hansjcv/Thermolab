function [Cname,phs_name] = exclude_phases(phs_name,Cname,Nsys,solution_model)
td    = init_thermo(phs_name,Cname,solution_model);
for ip = 1:length(phs_name)
    td(ip).dz(:) = 1/2;   
end
[~,Npc,pc_id] = tl_gibbs_energy(800+273,1e9,phs_name,td);
Cname  = Cname(Nsys>0);
for ip = 1:length(phs_name)
     exc_phs(ip) = sum(sum(Npc(Nsys==0,pc_id==ip)>0,1)==0)==0;
end
phs_name = phs_name(exc_phs==0);