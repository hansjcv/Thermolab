clear,addpath ../ ../Utilities/ ../Solutions/ ../EOS
T       = linspace(700,1000,4) + 273.15;
P       = 5e9;
R_cal   = 1.9858775;
R       = 8.3144;
Mw      = 18.01528e-3;
solids  = {'coe','diam','di','py','gr','jd','alm','O2'};
solvent = {'H2O'};
spcs    = {'H+','OH-','H2,aq',...
           'SiO2,aq','Si2O4,aq','Si3O6,aq','HSiO3-',...
           'CO3-2','CO2,aq','CO,aq','HCO3-','H2CO3,aq','METHANE,AQ','FORMATE,AQ','FORMIC-ACID,AQ','ACETATE,AQ','ACETATE,AQ','ACETIC-ACID,AQ',...
           'MgO,aq','Mg+2','MgOH+','Mg(HSiO3)+','Mg(SiO2)(HCO3)+','Mg(HCO3)+',...           
           'Ca+2','CaOH+','CaO,aq','Ca(HSiO3)+','Ca(HCO3)+','CaCO3,aq',...
           'AlO2-','AlO2(SiO2)-','HAlO2,aq','Al+3',...
           'Na+','NaOH,aq','NaHSiO3,aq','NaHCO3,aq','NaCO3-',...
           'Fe+2','Fe+3','Fe(CH3COO)+','FeOH+','FeO,aq','Fe(HCOO)+','Fe(CH3COO)2,aq','Fe(HSiO3)+'};
iCO2 = strcmp(spcs,'CO2,aq');
iCH4 = strcmp(spcs,'METHANE,AQ');
iH2 = strcmp(spcs,'H2,aq');
iO2 = strcmp(spcs,'O2,aq');
iCO = strcmp(spcs,'CO,aq');
iSi3O6 = strcmp(spcs,'Si3O6,aq');
iAlO2_SiO2m = strcmp(spcs,'AlO2(SiO2)-');
iNa = strcmp(spcs,'Na+');
%% Select dataset and solvent models
db_species    = 'DEW';
db_solids     = 'tc-ds633';
db_solvent    = 'tc-ds633';
solvent_model = 'ZD05';
dielect_model = 'S14';
gam_model     = 3;
%% Preprocess
phase   = [solids,solvent,spcs];
m_id    = [zeros(size(solids)),ones(size(solvent)),ones(size(spcs))*2];
for ip = 1:numel(phase)
    if m_id(ip) == 0
        phase{ip} = [phase{ip},',',db_solids];
    elseif m_id(ip) == 1
        phase{ip} = [phase{ip},',',db_solvent];
    elseif m_id(ip) == 2
        phase{ip} = [phase{ip},',',db_species];
    end    
end
%% Gibbs energy
[T2d,P2d]      = ndgrid(T,P);
[g0,Nphs]      = tl_gibbs_energy(T2d(:),P2d(:),phase);
[rho_w,eps_di] = water_props(T2d(:),P2d(:),solvent,solvent_model,dielect_model);
%% Reactions
niter          = 200;
max_Res        = 1e-14;
chrg           = Nphs(end,:);
v              = null(Nphs,'rational')';
no_dof         = sum(m_id==2)-size(v,1)-1;
disp(['Degrees of freedom = ' num2str(no_dof)]);
disp_reactions([solids,solvent,spcs],v');
act_solids = ones(numel(solids),length(T));
act_solids(end,:) = 10.^[-15.8,-13.6,-12,-10.4];
act_solids(strcmp(solids,'gr'),:) = [0.271,0.271,0.307,0.325];
act_solids(strcmp(solids,'py'),:) = [0.333,0.333,0.356,0.354];
act_solids(strcmp(solids,'alm'),:) =[ 0.396,0.396,0.336,0.321];
act_solids(strcmp(solids,'jd'),:) = [0.530,0.530,0.534,0.513];
act_solids(strcmp(solids,'di'),:) = [0.370,0.370,0.360,0.381];
logK = log10(exp(-(v(2,:)*g0).'/4.184./R_cal./T')); % To compare with DEW spreadsheet
%% Speciation calculations
for iT = 1:length(T)
    m0   = 0.9*ones(numel(phase(m_id==2)),1);
    m    = m0;
    gam  = ones(size(m));
    a_w  = 1;
    for iter = 1:niter
        act    = [act_solids(:,iT);a_w; m.*gam];        
        mu     = g0(:,iT)/8.3144/T(iT) + log(act);
        Res    = -[v*mu            ;chrg(m_id==2)*m];
        dRes   =  [v(:,m_id==2)./m';chrg(m_id==2)  ];        
        del_m  = dRes\Res;
        del_UR = 1./max(1,-del_m./(1/2.1*m));
        m      = m + del_UR.*del_m;
        z      = chrg(m_id==2)';    
        I      = m'*z.^2/2;        
        Lgam   = -log10(1+Mw*sum(m)); % large gamma (for 1 kg of water) conversion factor see Helgeson L gamma page 1293 eq 122, same as eq 2 in Miron 2016, assuming njw mole amount of water-solvent = 1
        a0     = 4.5e-8;
        A      = (1.82483e6)*sqrt(rho_w(iT)*1e-3) ./ (T(iT).*eps_di(iT)).^(3/2);
        B      = (50.2918649e8)*sqrt(rho_w(iT)*1e-3) ./ sqrt(T(iT).*eps_di(iT));
        b_gam  = 0.03*ones(size(z));
        % b_gam(z==0) = 0; % not better, but suggested in Huang
        b_gam(iCO2|iCH4|iH2|iO2|iCO) = -8.4495 + 0.01775*T(iT) - 7.5004e-6*(T(iT))^2;% gives bit closer fit
        % b_gam(iSi3O6) = 0.45;% my fit for Si
        % b_gam(iAlO2_SiO2m) = 0.3;
        % b_gam(iNa) = 0.001;
        if gam_model == 1 % Debye-Hückel
            lgam       = -(A*z.^2.*sqrt(I))./(1+a0*ones(size(z)).*B.*sqrt(I));
        elseif gam_model == 2 % Davies
            lgam = -A*z.^2.*(sqrt(I)./(1+sqrt(I))-0.2*I);
        elseif gam_model == 3 % HKF
            lgam    = -(A*z.^2.*sqrt(I))./(1+a0*ones(size(z)).*B.*sqrt(I)) + b_gam*I + Lgam; % probably molar scale see Helgeson L gamma page 1293 eq 122
        end
        gam = 10.^(lgam);
        % Activity of water
%         nu        = 2;
%         lam       = 1 + a0*B*sqrt(I); % lamda
%         sig       = 3./(a0^3*B.^3*sqrt(I.^3)).*(lam-1./lam-2*log(lam));  % sigma coefficient
%         phi       = -log(10)*sum(m(z~=0))./sum(m).*(A.*sqrt(I).*sig/3 + Lgam./(Mw*nu*I) - 0.03*I/2); % Osmotic coefficient
%         phi(I==0) = 0; % set osmotic coefficient to 0 for case of no electrolytes
%         a_w       = exp(-phi.*sum(m)*Mw); % Activity of water
        a_w    = (1/Mw)/(1/Mw+sum(m));
        if max(abs((m-m0)./m0).*100)<max_Res,break,end
        m0        = m;        
    end
    chk(iT)     = max(abs(Res(:)));
    pH(iT)      = -log10(m(strcmp(spcs,'H+')));    
    N_tot(iT,:) = Nphs(:,m_id==2)*m;   
    a_w_all(iT) = a_w;
    m_all(:,iT) = m;
end
%% Postprocess
[molality,id] = sort(m_all(:,2),'descend');
species  = spcs(id)';
table(species,molality)
disp(['pH = ' num2str(pH(2))])
disp(['aH2O = ' num2str(a_w_all(2))])
disp(['Degrees of freedom = ' num2str(no_dof)]);
disp(['chk = ',num2str(max(abs(chk(:))))])
Kessel_05_KFE_T_5GPa  = [700 800 900 1000];
Kessel_05_KFE_Si_5GPa = [0.57 2.29 5.87 11.23];
Kessel_05_KFE_Si_5GPa_err = [0.098 0.452 0.746 1.17];
Kessel_05_KFE_Al_5GPa = [0.09 0.4 0.43 1.21];
Kessel_05_KFE_Fe_5GPa = [0.05 0.1 0.09 0.21];
Kessel_05_KFE_Mg_5GPa = [0.09 0.24 0.16 0.35];
Kessel_05_KFE_Ca_5GPa = [0.19 0.12 0.50 1.04];
Kessel_05_KFE_Na_5GPa = [0.31 0.65 0.86 1.41];
figure(1),clf
subplot(321),plot(T-273.15,(N_tot(:,14)),'-'),axis([500 1100 0 13]),ylabel('molality of Si'),sgtitle('Eclogite at 5.0 GPa')
hold on
errorbar(Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Si_5GPa,Kessel_05_KFE_Si_5GPa_err,'o'),axis([500 1100 0 13]),ylabel('molality of Si'),sgtitle('Eclogite at 5.0 GPa')
subplot(322),plot(T-273.15,(N_tot(:,13)),'-',Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Al_5GPa,'o'),axis([500 1100 0 2]),ylabel('molality of Al'),legend('Thermolab','Kessel et al (2005) experiments','Location','NorthWest')
subplot(323),plot(T-273.15,(N_tot(:,26)),'-',Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Fe_5GPa,'o'),axis([500 1100 0 1]),ylabel('molality of Fe')
subplot(324),plot(T-273.15,(N_tot(:,12)),'-',Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Mg_5GPa,'o'),axis([500 1100 0 1]),ylabel('molality of Mg')
subplot(325),plot(T-273.15,(N_tot(:,20)),'-',Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Ca_5GPa,'o'),axis([500 1100 0 2]),ylabel('molality of Ca')
subplot(326),plot(T-273.15,(N_tot(:,11)),'-',Kessel_05_KFE_T_5GPa,Kessel_05_KFE_Na_5GPa,'o'),axis([500 1100 0 2]),ylabel('molality of Na')
% figure(2),clf
% iT = 3;
% yneg = [1 3 5 3 5 3 6 4 3 3];
% ypos = [2 5 3 5 2 5 2 2 5 5];
% xneg = [1 3 5 3 5 3 6 4 3 3];
% xpos = [2 5 3 5 2 5 2 2 5 5];
% errorbar(x,y,yneg,ypos,xneg,xpos,'o')
% plot([-1.5 1.5],[-1.5 1.5],...
%     log10(Kessel_05_KFE_Si_5GPa(iT)),log10(N_tot(iT,14)),'or',...
%     log10(Kessel_05_KFE_Al_5GPa(iT)),log10(N_tot(iT,13)),'or',...
%     log10(Kessel_05_KFE_Fe_5GPa(iT)),log10(N_tot(iT,26)),'or',...
%     log10(Kessel_05_KFE_Mg_5GPa(iT)),log10(N_tot(iT,12)),'or',...
%     log10(Kessel_05_KFE_Ca_5GPa(iT)),log10(N_tot(iT,20)),'or',...
%     log10(Kessel_05_KFE_Na_5GPa(iT)),log10(N_tot(iT,11)),'or'...
%     )
% axis([-1.5 1.5 -1.5 1.5]),axis square