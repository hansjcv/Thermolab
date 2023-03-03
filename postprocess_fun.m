function [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,g,alph] = postprocess_fun(T,P,td,alph,Npc,molm,p,pc_id,phases,solv_tol,feos,eps_model)
if ~exist('feos','var'),feos = 'CORK';end
if ~exist('eps_model','var'),eps_model = 'S14';end
% Numerics
delP     = 1e5;                                                          % for numerical differentiation
% Cluster analysis to find exsolved phases
[alph,Npc,p_out,pc_id] = cluster_p(alph,Npc,p,pc_id,solv_tol,phases);    % Compute exsolved, or numerically equivalent phases 
% For numerical derivatives
[rho_w_dP,eps_w_dP] = water_props(T,P+delP,phases,feos,eps_model);     % Calculate density of water
[rho_w,eps_w]       = water_props(T,P,phases,feos,eps_model);     % Calculate density of water
[g0,v0]       = tl_g0(T,P,td,rho_w,eps_w); % Calculate endmember Gibbs energies of each solution
[g0_dP,v0_dP] = tl_g0(T,P+delP,td,rho_w_dP,eps_w_dP);% Calculate endmember Gibbs energies of each solution at dP
g             = tl_gibbs_energy(T,P     ,phases,td,p_out,g0,v0,rho_w,eps_w);            % get Gibbs energy at P
g_P           = tl_gibbs_energy(T,P+delP,phases,td,p_out,g0_dP,v0_dP,rho_w_dP,eps_w_dP);% get Gibbs energy at P+dP
% Postprocessing
Vmol = (g_P-g)/delP;                                                     % Equation 44
Mmol = Npc'*molm;                                                        % Equation 45
rho  = Mmol./Vmol;                                                       % Equation 46
phim = alph/sum(alph);                                                   % phi mol
phi  = phim.*Vmol./(Vmol'*phim);                                         % Equation 47
phiw = phim.*Mmol./(Mmol'*phim);                                         % Equation 48
Cwt  = Npc.*repmat(molm,1,size(Npc,2))./repmat(Mmol',size(Npc,1),1);     % Equation 49
mu   = mu_fun(T,P,phases(pc_id),td(pc_id),p_out(pc_id),g0(pc_id),v0(pc_id),rho_w,eps_w);