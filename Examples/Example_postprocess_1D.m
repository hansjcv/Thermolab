clear,clf,addpath ../ ../Utilities/ ../Solutions/
runname = 'K32_2023_01_26_PTpath';
load(['linprog_run_' runname]);                                          % load linprog run data
solv_tol = 1;
fluid = 'H2O,tc-ds55';
grt_frac = [];
grt_comp = [];
molm  = molmass_fun(Cname);
phs_modes = zeros(length(T(:)),length(phs_name));
cnt = 0;
for iPT = 1:length(T(:))    
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,g,alph] = postprocess_fun(T(iPT),P(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,fluid,solv_tol);         
    phs_modes(iPT,pc_id) = phi;
    if sum(pc_id==3)>0
        phi_grt(iPT) = phi(pc_id==3);
        Ngrt(:,iPT)  = Npc(:,pc_id==3);
        grt_frac     = [grt_frac  phi_grt(iPT)];
        grt_comp     = [grt_comp Ngrt(:,iPT)];
    else
        phi_grt(iPT) = 0;
        Ngrt(:,iPT)  = zeros(length(Cname),1);
    end
end
vol_grt = cumsum(phi_grt);
solid_id    = find(~strcmp(phs_name,fluid));
solid_names = phs_name(solid_id);
vol_frac_solids = phs_modes(:,solid_id)./sum(phs_modes(:,solid_id),2);
solid_names = solid_names(sum(vol_frac_solids)>0);
vol_frac_solids = vol_frac_solids(:,sum(vol_frac_solids)>0);
vol_frac_solids(:,strcmp(solid_names,'Garnet')) = vol_grt;
vol_frac_solids = vol_frac_solids./sum(vol_frac_solids,2);
figure(1),
plot(T-273.15,cumsum(phi_grt),'o-');
figure(2),
plot(P/1e9,Ngrt(strcmp(Cname,'Fe'),:),'o-',P/1e9,Ngrt(strcmp(Cname,'Mg'),:),'o-',P/1e9,Ngrt(strcmp(Cname,'Ca'),:),'o-')
figure(3)
plot(cumsum(grt_frac),grt_comp([4,5,6,7],:),'o-');legend(Cname([4,5,6,7]))
figure(4),colormap jet
area(T-273.15,vol_frac_solids,'FaceColor','flat'),axis tight,legend(solid_names)