clear,addpath ../../ ../../EOS ../../Solutions/ ../../Utilities
R = 8.3145;
proportions = [];
sitefractions = [];
mu0           = [];
GibbsEnergies = [];
g_ideal = [];
g_non_ideal = [];
Minerals = [];
for case_id = 1:13
%% Load data
if 1 == case_id
    amp_data;
    phase = {'Amphibole'};
elseif 2 == case_id  
    cpx_data;
    phase = {'Clinopyroxene'};
elseif 3 == case_id
    grt_data;
    phase = {'Garnet'};
elseif 4 == case_id
    ep_data;
    phase = {'Epidote'};
elseif 5 == case_id
    mu_data;
    phase = {'Muscovite'};
elseif 6 == case_id
    ol_data;
    phase = {'Olivine'};
elseif 7 == case_id
    bt_data;
    phase = {'Biotite'};
elseif 8 == case_id
    pl4tr_data;
    phase = {'Feldspar'};
elseif 9 == case_id
    opx_data;
    phase = {'Orthopyroxene'};
elseif 10 == case_id
    ilm_data;
    phase = {'Ilmenite'};
elseif 11 == case_id
    spl_data;
    phase = {'Spinel'};
elseif 12 == case_id
    liq_data;
    phase = {'Melt(G25)'};
elseif 13 == case_id
    fl_data;
    phase = {'Fluid(G25)'};
end
chk_p = [];
chk_z = [];
chk_g0 = [];
chk_a_id = [];
chk_gid = [];
chk_gam = [];
chk_g = [];
chk_RTlngam = [];
chk_RTlna_id = [];
chk_RTlna = [];
%% Calculate Gibbs and compare
for i = 1:length(T_tc)
    Temp  = T_tc(i);
    Pres  = P_tc(i);
    Cname = {'Si','Al','Ca','Mg','Fe','K','Na','Ti','H','Cr','O'};
    td    = init_thermo(phase,Cname,'Igneous');
    [T,P] = ndgrid(Temp,Pres);
    [rho_w,eps_di]  = water_props(T(:),P(:),phase,'PS94','S14');
    [g0,v0] = tl_g0(T(:),P(:),td,rho_w,eps_di,1); % use option 1 to use 'apparent Gibbs energy'
    p{1} = p_tc(i,:);
    [g,Npc,pc_id,p1,z,g_id,g_nid] = tl_gibbs_energy(T(:),P(:),phase,td,p,g0,v0);
    [mu,a,RTlngam] = tl_chemical_potential(T(:),P(:),td,p{1},g0);
    z_indep         = z_tc(i,td.site_var);
    p2              = (td.p_from_z_cons*[ones(1,size(z_indep,1)); z_indep'])';
    chk_p(i)        = max(abs(p2-p_tc(i,:)));
    chk_z(i)        = max(abs(z-z_tc(i,:)));
    chk_g0(i)       = max(abs(g0{1}-mu0_tc(i,:)*1e3));
    chk_a_id(i)     = max(abs(a-a_id_tc(i,:)));
    chk_gid(i)      = max(abs(g_id-R*T*log(a_id_tc(i,:) + double(a_id_tc(i,:)==0))*p_tc(i,:)'));
    chk_gnid(i)     = max(abs(g_nid-R*Temp*log(gam_tc(i,:))*p_tc(i,:)'));
    chk_gam(i)      = max(abs(gam_tc(i,:) - exp(RTlngam/Temp/R)));
    chk_g(i)        = max(abs(g-g_tc(i)*1e3));
    chk_RTlngam(i)  = max(abs(R*Temp*log(gam_tc(i,:)) - (RTlngam)));
    chk_RTlna_id(i) = max(abs(R*Temp*log(a_id_tc(i,:)) - (R*Temp*log(a))));
    chk_RTlna(i)    = max(abs(RTlna_tc(i,:) - (R*Temp*log(a)+RTlngam)));
end
proportions = [proportions;max(chk_p)'];
sitefractions = [sitefractions;max(chk_z)'];
mu0 = [mu0;max(chk_g0)'];
g_ideal = [g_ideal;max(chk_gid)];
g_non_ideal = [g_non_ideal;max(chk_gnid)];
GibbsEnergies = [GibbsEnergies;max(chk_g)'];
Minerals = [Minerals;phase];
end
%% Show results
disp('This table shows maximum difference in Gibbs energies in Joules')
table(GibbsEnergies,mu0,g_ideal,g_non_ideal,proportions,sitefractions,'RowNames',Minerals)