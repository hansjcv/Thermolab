clear,addpath ../../ ../../Solutions ../../EOS ../../Utilities/
R = 8.3145;
proportions = [];
sitefractions = [];
mu0           = [];
GibbsEnergies = [];
g_ideal = [];
g_non_ideal = [];
Minerals = [];
for case_id = 12
    %% Load data
    if 1 == case_id
        atg_data;
        phase = {'Antigorite'};
    elseif 2 == case_id
        chl_data;
        phase = {'Chlorite'};
    elseif 3 == case_id
        br_data;
        phase = {'Brucite'};
    elseif 4 == case_id
        ol_data;
        phase = {'Olivine'};
    elseif 5 == case_id
        ch_data;
        phase = {'Clinohumite'};
    elseif 6 == case_id
        po_data;
        phase = {'Pyrrhotite'};
    elseif 7 == case_id
        spi_data;
        phase = {'Spinel'};
    elseif 8 == case_id
        opx_data;
        phase = {'Orthopyroxene'};
    elseif 9 == case_id
        ta_data;
        phase = {'Talc'};
    elseif 10 == case_id
        g_data;
        phase = {'Garnet'};
    elseif 11 == case_id
        anth_data;
        phase = {'Anthophyllite'};
    elseif 12 == case_id
        fl_data;
        phase = {'Fluid'};
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
        Cname = {'Si','Ti','Al','Fe','Mn','Mg','Ca','Na','K','S','H','O'};
        td    = init_thermo(phase,Cname,'Ultramafic');
        [T,P] = ndgrid(Temp,Pres);
        [rho_w,eps_di]  = water_props(T(:),P(:),phase,'PS94','S14');
        [g0,v0] = tl_g0(T(:),P(:),td,rho_w,eps_di,1); % use option 1 to use 'apparent Gibbs energy'
        p{1} = p_tc(i,:);
        [g,Npc,pc_id,p1,z,g_id,g_nid] = tl_gibbs_energy(T(:),P(:),phase,td,p,g0,v0);
        [mu,a,RTlngam] = tl_chemical_potential(T(:),P(:),td,p{1},g0);
        z_indep      = z(td.site_var);%z_tc(i,td.site_var);
        p2           = (td.p_from_z_cons*[ones(1,size(z_indep,1)); z_indep'])';       
        chk_z(i)     = max(abs(z-z_tc(i,:)));
        chk_p(i)     = max(abs(p2-p_tc(i,:)));        
        chk_g0(i)    = max(abs(g0{1}-mu0_tc(i,:)*1e3));
        chk_a_id(i)  = max(abs(a-a_id_tc(i,:)));
        chk_gid(i)   = max(abs(g_id-R*T*log(a_id_tc(i,:) + double(a_id_tc(i,:)==0))*p_tc(i,:)'));
        chk_gnid(i)   = max(abs(g_nid-R*Temp*log(gam_tc(i,:)+ double(gam_tc(i,:)==0))*p_tc(i,:)'));
        chk_gam(i)   = max(abs(gam_tc(i,:) - exp(RTlngam/Temp/R)));
        chk_g(i)     = max(abs(g-g_tc(i)*1e3));
        chk_RTlngam(i) = max(abs(R*Temp*log(gam_tc(i,:)) - (RTlngam)));
        chk_RTlna_id(i) = max(abs(R*Temp*log(a_id_tc(i,:)) - (R*Temp*log(a))));
        chk_RTlna(i) = max(abs(RTlna_tc(i,:) - (R*Temp*log(a)+RTlngam)));
    end
    proportions   = [proportions;max(chk_p)'];
    sitefractions = [sitefractions;max(chk_z)'];
    mu0           = [mu0;max(chk_g0)'];
    g_ideal       = [g_ideal;max(chk_gid)];
    g_non_ideal   = [g_non_ideal;max(chk_gnid)];   
    GibbsEnergies = [GibbsEnergies;max(chk_g)'];
    Minerals      = [Minerals;phase];
end
%% Show results
disp('This table shows maximum difference in Gibbs energies in Joules comparing THERMOCALC output with Thermolab')
table(GibbsEnergies,mu0,g_ideal,g_non_ideal,proportions,sitefractions,'RowNames',Minerals)
