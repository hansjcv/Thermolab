clear,addpath ../ ../Solutions ../Utilities/ ../EOS/
%% Input
T       = 995.3340;                                                                 % T in Kelvin
P       = 0.2e9;                                                                    % P in Pascal
solname = 'Amphibole';                                                              % Name of the phase
solfile = 'solution_models_H18';                                                    % Solution model file name
Cname   = {'Si'  ,'Ti'  ,'Al'  ,'Ca',  'Fe',  'Mg',  'Na'  ,'K'   ,'H'   ,'O'    }; % System component names
Cmol    = [7.0006,0.1054,1.3603,1.6786,1.8705,2.9542,0.3382,0.0237,1.8946,24.0500]; % Composition in mol
%% Load thermodynamic data
[p_id,st,site_id,mtpl,mod_id,alp,w,n_em,chg,em_data,dGex,EOS,CEOS,mcoef,Gr,p_name] = init_phase(solfile,solname,Cname);
%% Water properties
rho_w   = rho_H2O(T,P,'ZD05');       % Water density (ZD05= Zhang & Duan (2005))
eps_di  = eps_H2O(T,P,rho_w,'S14');  % Dielectric constant (S14 = Sverjensky et al.2014)
%% Endmember Gibbs energies
for k = 1:length(em_data)
    Gjp = zeros(1,size(em_data{k},1));                                     % Initialize G of make endmembers
    for jp = 1:size(em_data{k},1)                                          % For endmembers made of multiple phases
        [SdT,S]   = intSdT(T,P,em_data{k}(jp,:),CEOS{k}(jp));              % SdT integral, (eq. 10&11)
        [VdP,V]   = intVdP(T,P,em_data{k}(jp,:), EOS{k}(jp));              % VdP integral, (eq. 12)
        Gexc      = g0_exc(T,P,em_data{k}(jp,:),dGex{k}(jp),rho_w,eps_di); % dG of any fitting excess energy (for eq. 12)
        Gjp(:,jp) = Gr{k}(jp) - SdT + VdP + Gexc;                          % Eq. 12
    end
    g0(k) = mcoef{k}*Gjp';                                                 % Eq. 17
end
%% Prepare conversion matrices 
zt                               = st2zt(st,site_id);               % Site fraction table, (eq. 25)
[p_from_c_cons,icomp_indep,isite_indep] = comp2prop(n_em',zt');     % C variables, p to C conversion (Eq.41-43 & App. A)
%% Grid for order-disorder
[z_od1,z_od2]                     = ndgrid(linspace(0,1,100),linspace(0,1,100));    % all possible order-disorder states
z_od                              = [z_od1(:),z_od2(:)];                            % order-disorder state matrix
%% Proportions and site fractions
p  = (p_from_c_cons*[ones(length(z_od),1) repmat(Cmol(icomp_indep),length(z_od),1) z_od]')'; % Eq. 43 & App.A
z  = p*zt;                                                                                   % Eq. 22 (site fractions)
%% Numerical
z_tol = 1e-10;z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;   % round off errors
badz  = min(z,[],2)<0|max(z,[],2)>1;p(badz,:) = [];z(badz,:) = [];             % remove unrealistic compositions
%% Gibbs energy of mixing
g_mech = p*g0';                                                                             % Mechanical mixing, eq. 19
g_id   = (T*8.3145*(mtpl*(z'.*log(z'+double(z'==0))- (zt'.*log(zt'+double(zt'==0)))*p')))'; % Ideal mixing, eq. 20
g_nid  = gnid(T,P,p,mod_id,V,alp,w,rho_w,eps_di,chg);                                       % Non-ideal mixing. eq. 21
g      = g_mech + g_id + g_nid;                                                             % Gibbs energy in Joule/mol, eq. 18
%% Find stable configuration (homogeneous equilibrium)
[gmin,id] = min(g);            % minimum of Gibbs energy
p(id,:)                        % corresponding endmember proportions
%% Show order-disorder reactions
disp_reactions(p_name,null(n_em','r'))