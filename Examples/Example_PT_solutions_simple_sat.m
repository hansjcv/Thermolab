clear,clf,addpath ../ ../EOS ../Solutions ../Utilities
runname  = 'test';
lp_run   = 1;
T        = linspace(300,800,55) + 273.15;
P        = linspace(0.1,3.5,56)*1e9;
Cname    = {'SiO2','FeO','MgO' ,'H2O'};
molm     = molmass_oxy(Cname);
Nsys     = [34   , 3     45  , 31];
db_name  = 'tc-ds55';
fluid    = 'H2O,tc-ds55';
ndim     = 2;
ps_type  = 'P-T';
solv_tol2 = 1;
sol_name = {'Antigorite','Olivine','Brucite','Talc','Orthopyroxene'};
pur_name = {'anth,tc-ds55'};%get_pure_phases(db_name,Cname,is_liq,is_gas,is_aq);%{'q,tc-ds55','an,tc-ds55','gr,tc-ds55','ky,tc-ds55'};
phs_name = [sol_name,fluid,pur_name];
td       = init_thermo(phs_name,Cname,'solution_models_EF21');
p        = props_generate(td);     % generate endmember proportions
[T2d,P2d] = ndgrid(T,P);
[g,Npc,pc_id] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td,p);
g = g-g(pc_id==find(strcmp(phs_name,'H2O,tc-ds55')),:).*Npc(strcmp(Cname,'H2O'),:)';
LB = zeros(1,size(g,1));
if lp_run == 1
    % Minimization refinement
    tic
    for iPT = 1:length(T2d(:))
        [alph(iPT,:),fval,exitflag,output,lambda] = linprog(g(:,iPT),[],[],Npc(1:end-1,:),Nsys(1:end-1),LB);        
        [pc_id_out,phi_out,Cwt_out,Npc_out,rho_out,mu_out,p_out,phiw_out,g_out,alph_out] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph(iPT,:)',Npc,molm,p,pc_id,phs_name,solv_tol2,'CORK','S14');
        asm_id(iPT,1:length(unique(pc_id(alph(iPT,:)>0)))) = unique(pc_id(alph(iPT,:)>0));        
        disp(iPT/length(T2d(:)))
    end
    toc
    save(['linprog_run_' runname],'-v7.3');
else
    load(['linprog_run_' runname]);
end
figure(1),tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name,0,[0 0],10);