clear,figure(2),clf,colormap jet,addpath ../
run_name = 'soapstone_2022_02_14_n1500_nz15_full';
load(['linprog_run_' run_name]);
options.eps_solvus = 10;
delP = 1e5;
delc = 1e-3;
ET   = struct([]);
melt = 'Melt(W07)';
phi   = zeros(length(T2d(:)),length(phs_name));
molm = molmass_fun(Cname);
% Step 4:Postprocessing
for iPT = 1:length(T2d(:))
    [N_mol,Npc,p,pc_id] = cluster_p(alph_all{iPT},Npc_all{iPT},p_all{iPT},pc_id_all{iPT},options.eps_solvus,phs_name);
    [g0_p,V0_p]         = tl_g0(T2d(iPT),P2d(iPT)+delP,td,rho_w,eps_di);
    [g0  ,V0  ]         = tl_g0(T2d(iPT),P2d(iPT)     ,td,rho_w,eps_di);
    g                   = tl_gibbs_energy(T2d(iPT),P2d(iPT),phs_name(unique(pc_id)),td(unique(pc_id)),p(unique(pc_id)),g0(unique(pc_id)),V0(unique(pc_id)),rho_w,eps_di);
    g_dP                = tl_gibbs_energy(T2d(iPT),P2d(iPT)+delP,phs_name(unique(pc_id)),td(unique(pc_id)),p(unique(pc_id)),g0_p(unique(pc_id)),V0_p(unique(pc_id)),rho_w,eps_di);
    V_mol               = (g_dP-g)/delP;
    M_mol               = molm'*Npc; % calculate molar masses of the phases
    phi_mol             = N_mol/sum(N_mol); % calculate mol fraction
    ET(iPT).Npc         = Npc;
    ET(iPT).pc_id       = pc_id;    
    ET(iPT).phi         = phi_mol.*V_mol/(V_mol'*phi_mol); % volume fraction    
    ET(iPT).rho         = M_mol'./V_mol;% kg/m3 calculate density of all phases        
    ET(iPT).C_wt        = ET(iPT).Npc.*molm./(molm'*ET(iPT).Npc);% weight composition        
    ET(iPT).phi_wt      = phi_mol.*M_mol'/sum(phi_mol.*M_mol');% weight fraction    
    fluid_id            =   strcmp(phs_name(ET(iPT).pc_id),fluid)|strcmp(phs_name(ET(iPT).pc_id),melt);
    solid_id            = ~(strcmp(phs_name(ET(iPT).pc_id),fluid)|strcmp(phs_name(ET(iPT).pc_id),melt));    
    p_fluid{1}          = [p{fluid_id}(1); 1-p{fluid_id}(1)]';                                   % proportions of endmembers Ab-Or
    p_fluid_dC{1}       = [p{fluid_id}(1)+delc; 1-(p{fluid_id}(1)+delc)]';                       % proportions of endmembers Ab-Or plus delta C
    g_fl                = tl_gibbs_energy(T,P,phs_name(fluid_id),td(fluid_id),p_fluid,g0(fluid_id),V0(fluid_id)); % Gibbs energy + delta G, for numerical derivative    
    g_fl_dC             = tl_gibbs_energy(T,P,phs_name(fluid_id),td(fluid_id),p_fluid_dC,g0(fluid_id),V0(fluid_id)); % Gibbs energy + delta G, for numerical derivative    
    mu_fluid(iPT)       = (g_fl_dC-g_fl)'/delc;                                                  % mu1 - mu2; numerical derivative dg/dc
    mu_all(:,iPT)       = Npc'\g;
    mu_fluid2(iPT)      = Npc(end-3,Npc(end-3,:)>0)'\g(Npc(end-3,:)>0);    
    rhos_tab(iPT)       = ET(iPT).rho(solid_id)'*ET(iPT).phi(solid_id)/sum(ET(iPT).phi(solid_id));
    rhof_tab(iPT)       = ET(iPT).rho(fluid_id)'*ET(iPT).phi(fluid_id)/sum(ET(iPT).phi(fluid_id));    
    cwt_solid(:,iPT)    = ET(iPT).C_wt(:,solid_id)*ET(iPT).phi_wt(solid_id)/sum(ET(iPT).C_wt(:,solid_id)*ET(iPT).phi_wt(solid_id));
    cwt_fluid(:,iPT)    = ET(iPT).C_wt(:,fluid_id)*ET(iPT).phi_wt(fluid_id)/sum(ET(iPT).C_wt(:,fluid_id)*ET(iPT).phi_wt(fluid_id));    
    Cs_tab(iPT)         = cwt_solid(strcmp(Cname,'C'),iPT);
    Mg_tab(iPT)         = cwt_solid(strcmp(Cname,'Mg'),iPT);        
    Cf_tab(iPT)         = cwt_fluid(strcmp(Cname,'C'),iPT);
    phi(iPT,pc_id)      = ET(iPT).phi;
    N_mol_all{iPT}      = N_mol;
    disp(iPT/length(T2d(:)))
end
mu_tab = mu_fluid;
solid_id    = find(~strcmp(phs_name,fluid));
solid_names = phs_name(solid_id);
vol_frac_solids = phi(:,solid_id)./sum(phi(:,solid_id),2);
solid_names = solid_names(sum(vol_frac_solids)>0);
[solid_names,abc_id] = sort(solid_names);
vol_frac_solids = vol_frac_solids(:,sum(vol_frac_solids)>0);
vol_frac_solids = vol_frac_solids(:,abc_id);
subplot(211),area(Cs_tab,vol_frac_solids),legend(solid_names),axis tight,axis square
set(gca,'FontSize',14)
subplot(212),plot(Cs_tab,Cf_tab,'*-','LineWidth',0.5),axis tight,axis square
xlabel('C solid')
ylabel('C fluid')
set(gca,'FontSize',14)
save(['lookup_' run_name], 'Cf_tab', 'Cs_tab', 'vol_frac_solids', 'solid_names', 'rhos_tab', 'rhof_tab', 'Mg_tab','mu_tab')