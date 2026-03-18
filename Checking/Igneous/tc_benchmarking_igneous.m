clear,addpath ../../ ../../EOS ../../Solutions/ ../../Utilities
for i_case = 1:13
case_id = i_case;
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

for i = 1:length(T_tc)
    Temp  = T_tc(i);
    Pres  = P_tc(i);
    Cname = {'Si','Al','Ca','Mg','Fe','K','Na','Ti','H','Cr','O'};
    td    = init_thermo(phase,Cname,'Igneous');
    [T,P] = ndgrid(Temp,Pres);
    [rho_w,eps_di]  = water_props(T(:),P(:),phase,'PS94','S14');
    [g0,v0] = tl_g0(T(:),P(:),td,rho_w,eps_di,1); % use option 1 to use 'apparent Gibbs energy'
    p{1} = p_tc(i,:);%[-0.21756,0.35166,0.56365,-0.038051,0.0037738,0.30484,-0.28228,-0.00075077,0.21444,0.076988,0.023286];        
    [g,Npc,pc_id,p1,z,g_id,g_nid] = tl_gibbs_energy(T(:),P(:),phase,td,p,g0,v0);    
    % p2 = (td.p_from_z_cons*[ones(1,size(td.site_var,1)); z(td.site_var)'])';
    [mu,a,RTlngam] = tl_chemical_potential(T(:),P(:),phase,td,p{1},Cname,g0,v0);
    chk_z(i)     = max(abs(z-z_tc(i,:)));
    chk_g0(i)    = max(abs(g0{1}-mu0_tc(i,:)*1e3));
    chk_gid(i)   = max(abs(g_id-8.3144*T*log(a_id_tc(i,:))*p_tc(i,:)'));
    chk_a_id(i)  = max(abs(a-a_id_tc(i,:)));
    chk_gam(i)   = max(abs(gam_tc(i,:) - exp(RTlngam/Temp/8.3144)));
    chk_g(i)     = max(abs(g-g_tc(i)*1e3));
    chk_RTlna(i) = max(abs(RTlna_tc(i,:)*1e3 - (8.3144*Temp*log(a)+RTlngam)));
end

chk_z
chk_g0
chk_a_id
chk_gam
chk_RTlna
chk_gid
chk_g 
end