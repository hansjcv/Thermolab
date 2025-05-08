clear,clf, addpath ../ ../Utilities/ ../EOS
T = linspace(200,900,100) + 273.15;
P = 2.0e8;
[T2d,P2d] = ndgrid(T,P);
for i_dataset = 1:2
    if i_dataset == 1
        spc       = {'HF,aq,supcrt','mu,tc-ds55','FLUORMUSCOVITE,ZS91','H2O,tc-ds55'};
        spc_buff  = {'an,tc-ds55','HF,aq,supcrt','FLUORITE,supcrt','sill,tc-ds55','q,tc-ds55','H2O,tc-ds55'};   
        % spc       = {'HF,aq,supcrt','MUSCOVITE,ZS91','FLUORMUSCOVITE,ZS91','H2O,tc-ds55'};
    elseif i_dataset == 2
        spc       = {'HF,aq,supcrt','mu,tc-ds55','fmu,ZS91','H2O,tc-ds55'};
        spc_buff  = {'ANORTHITE,ZS91','HF,aq,supcrt','FLUORITE,supcrt','SILLIMANITE,ZS91','A-QUARTZ,ZS91','H2O,tc-ds55'};   
    end
    
    
    Cname     = {'Si','Ca','K','F','Al','H','O','e'};
    %% Thermolab Gibbs energies
    td        = init_thermo(spc,Cname);
    p         = props_generate(td);
    [rho_w,eps_di] = water_props(T2d(:),P2d(:),spc,'IAPWS','S14');
    [g0,v0]        = tl_g0(T2d(:),P2d(:),td,rho_w,eps_di);
    [g,Npc]        = tl_gibbs_energy(T2d(:),P2d(:),spc,td,p,g0,v0);
    %% Thermolab Gibbs energies buffer
    td_buff        = init_thermo(spc_buff,Cname);
    p_buff         = props_generate(td_buff);
    [rho_w_buff,eps_di_buff] = water_props(T2d(:),P2d(:),spc_buff,'IAPWS','S14');
    [g0_buff,v0_buff]        = tl_g0(T2d(:),P2d(:),td_buff,rho_w_buff,eps_di_buff);
    [g_buff,Npc_buff]        = tl_gibbs_energy(T2d(:),P2d(:),spc_buff,td_buff,p_buff,g0_buff,v0_buff);
    g_mu        = tl_gibbs_energy(T2d(:),P2d(:),{'mu,tc-ds55','MUSCOVITE,ZS91'});
    g_muf        = tl_gibbs_energy(T2d(:),P2d(:),{'FLUORMUSCOVITE,ZS91','fmu,ZS91'});
    %% Reactions
    react_buff = null(Npc_buff,'r');
    react = null(Npc,'r');
    disp_reactions(spc,react);
    disp_reactions(spc_buff,react_buff);
    dg_buff  = react_buff'*g_buff;
    Keq_buff = exp(-dg_buff./8.3144./T);
    logaHF = log10(Keq_buff)/2;
    aHF    = 10.^(logaHF);
    g(1,:) = g(1,:) + 8.3144*T.*log(aHF);
    dg  = react'*g;
    Keq = exp(-dg./8.3144./T);
    %% Plotting
    subplot(211),hold on,plot(T-273.15,log10(Keq))    
    leg_names{i_dataset} = disp_reactions(spc,react);
    subplot(212),hold on,plot(T-273.15,log10(Keq_buff))
    leg_names2{i_dataset} = disp_reactions(spc_buff,react_buff);
    hold on
end
T_ZS91    = [462 462 500 590 612 612 640 690];
logK_ZS91 = [1.24 2.36 1.44 0.89 0.65 0.81 0.48 0.46];
logKbuff_ZS91 = [3.8 3.8 3.26 2.15 1.91 1.91 1.63 1.16];
subplot(211),hold on,plot(T_ZS91,logK_ZS91,'o'),legend([leg_names,'Zhu and Sverjensky (1991)']),title('Comparison fluormuscovite endmember')
axis([200 900 -1 3]),xlabel('T(\circC)'),ylabel('log K')
subplot(212),plot(T_ZS91,logKbuff_ZS91,'o'),xlabel('T(\circC)'),ylabel('log K'),legend([leg_names2,'Zhu and Sverjensky (1991)'])