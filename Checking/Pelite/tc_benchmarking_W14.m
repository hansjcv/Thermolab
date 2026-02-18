clear,addpath ../../ ../../Solutions ../../EOS ../../Utilities/
case_id = 17;
if 1 == case_id
    st_data;
    phase = {'Staurolite'};
elseif 2 == case_id  
    ma_data;
    phase = {'Margarite'};
elseif 3 == case_id
    g_data;
    phase = {'Garnet'};
elseif 4 == case_id
    ep_data;
    phase = {'Epidote'};
elseif 5 == case_id
    mu_data;
    phase = {'Muscovite'};
elseif 6 == case_id
    cd_data;
    phase = {'Cordierite'};
elseif 7 == case_id
    bi_data;
    phase = {'Biotite'};
elseif 8 == case_id
    pl4tr_data;
    phase = {'Feldspar'};
elseif 9 == case_id
    chl_data;
    phase = {'Chlorite'};
elseif 10 == case_id
    opx_data;
    phase = {'Orthopyroxene'};
elseif 11 == case_id
    ilm_data;
    phase = {'Ilmenite'};
elseif 12 == case_id
    mt1_data;
    phase = {'Magnetite'};
elseif 13 == case_id
    ctd_data;
    phase = {'Chloritoid'};      
elseif 14 == case_id
    sp_data;
    phase = {'Spinel'};     
elseif 15 == case_id
    sa_data;
    phase = {'Sapphirine'};     
elseif 16 == case_id
    ilmm_data;
    phase = {'Ilmenite-Mn'};   
elseif 17 == case_id
    liq_data;
    phase = {'Melt(W14)'};    
end

for i = 1:length(T_tc)
    Temp  = T_tc(i);
    Pres  = P_tc(i);
    Cname = {'Si','Ti','Al','Fe','Mn','Mg','Ca','Na','K','H','O'};
    td    = init_thermo(phase,Cname,'Metapelite');
    [T,P] = ndgrid(Temp,Pres);
    [rho_w,eps_di]  = water_props(T(:),P(:),phase,'PS94','S14');
    [g0,v0] = tl_g0(T(:),P(:),td,rho_w,eps_di,1); % use option 1 to use 'apparent Gibbs energy'
    p{1} = p_tc(i,:);%[-0.21756,0.35166,0.56365,-0.038051,0.0037738,0.30484,-0.28228,-0.00075077,0.21444,0.076988,0.023286];
    [g,Npc,pc_id,p1,z,g_id,g_nid] = tl_gibbs_energy(T(:),P(:),phase,td,p,g0,v0);
    [mu,a,RTlngam] = tl_chemical_potential(T(:),P(:),phase,td,p{1},Cname,g0,v0);
    z_indep      = z_tc(i,td.site_var);
    p2           = (td.p_from_z_cons*[ones(1,size(z_indep,1)); z_indep'])';
    chk_p(i)     = max(abs(p2-p_tc(i,:)));
    chk_z(i)     = max(abs(z-z_tc(i,:)));
    chk_g0(i)    = max(abs(g0{1}-mu0_tc(i,:)*1e3));
    chk_a_id(i)  = max(abs(a-a_id_tc(i,:)));
    chk_gam(i)   = max(abs(gam_tc(i,:) - exp(RTlngam/Temp/8.3145)));
    chk_g(i)     = max(abs(g-g_tc(i)*1e3));
    chk_RTlngam(i) = max(abs(8.3145*Temp*log(gam_tc(i,:)) - (RTlngam)));
    chk_RTlna_id(i) = max(abs(8.3145*Temp*log(a_id_tc(i,:)) - (8.3145*Temp*log(a))));
    chk_RTlna(i) = max(abs(RTlna_tc(i,:) - (8.3145*Temp*log(a)+RTlngam)));
end

chk_p
chk_z
chk_g0
chk_a_id
chk_gam
chk_RTlngam
chk_RTlna_id
chk_RTlna
chk_g


%% Notes

% check gamma for biotite 'tbi,tc-ds62' - difference on 3rd decimal place

