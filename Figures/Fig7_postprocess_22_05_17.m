clear,addpath ../
runname = 'example_soap';
load([runname '_linprog']);                                              % load linprog run data
molm = molmass_fun(Cname);                                               % get molar mass of the components
% Numerics
delP     = 1e5;                                                          % for numerical differentiation
delc     = 1e-5;                                                         % for numerical differentiation
solv_tol = 2;                                                            % Tolerance for recognizing exsolved phases
% Cluster analysis to find exsolved phases
[alph,Nphs,p_out,pc_id] = cluster_p(alph,Npc,p,pc_id,solv_tol,phases);    % Compute exsolved, or numerically equivalent phases 
% For numerical derivatives
[rho_w   ,eps_w   ] = water_props(T,P     ,phases,'ZD05','S14');          % Calculate density of water
[rho_w_dP,eps_w_dP] = water_props(T,P+delP,phases,'ZD05','S14');          % Calculate density of water
[g0,v0]       = tl_g0(T,P,td,rho_w,eps_w);            % Calculate endmember Gibbs energies of each solution
[g0_dP,v0_dP] = tl_g0(T,P+delP,td,rho_w_dP,eps_w_dP); % Calculate endmember Gibbs energies of each solution at dP
g             = tl_gibbs_energy(T,P     ,phases,td,p_out,g0,v0,rho_w,eps_w);            % get Gibbs energy at P
g_P           = tl_gibbs_energy(T,P+delP,phases,td,p_out,g0_dP,v0_dP,rho_w_dP,eps_w_dP);% get Gibbs energy at P+dP
% Postprocessing
Vmol     = (g_P-g)/delP;                                                     % Equation 54
Mmol     = Nphs'*molm;                                                       % Equation 56
rho      = Mmol./Vmol;                                                       % Equation 57
phim     = alph/sum(alph);                                                   % Equation 52
phi      = phim.*Vmol./(Vmol'*phim);                                         % Equation 53
phiw     = phim.*Mmol./(Mmol'*phim);                                         % Equation 55
Cwt      = Nphs.*repmat(molm,1,size(Nphs,2))./repmat(Mmol',size(Nphs,1),1);  % Equation 58
fluid_id =  strcmp(phases(pc_id),fluid);                                     % Find index of fluid
solid_id = ~fluid_id;                                                        % Find index of solids
rhos     = rho(solid_id)'*phi(solid_id)/sum(phi(solid_id));                  % Equation 59
rhof     = rho(fluid_id)'*phi(fluid_id)/sum(phi(fluid_id));                  % Fluid density
% Chemical potentials
for iphase = 1:numel(g)
    Gphase  = g(iphase);
    mu_last = Gphase;
    nc      = sum(p_out{pc_id(iphase)}>0);
    if nc > 1         
        for ic = 1:nc-1
            c_ind                        = find(p_out{pc_id(iphase)}>0);
            p1                           = p_out;
            p1{pc_id(iphase)}(c_ind(ic)) = p_out{pc_id(iphase)}(c_ind(ic)) + delc;
            p1{pc_id(iphase)}(c_ind(nc)) = p_out{pc_id(iphase)}(c_ind(nc)) - delc;
            g_dc                         = tl_gibbs_energy(T,P,phases,td,p1,g0,v0,rho_w,eps_w);
            Gphase1                      = g_dc(iphase);
            dGdc1                        = (Gphase1-Gphase)/delc;
            mu_last                      = mu_last - p_out{pc_id(iphase)}(c_ind(ic)).*dGdc1; % Equation 60
            mu{iphase}(ic)               = dGdc1;
        end
        for ic = 1:nc-1
            mu{iphase}(ic) = mu{iphase}(ic) + mu_last;                                       % Equation 61
        end
    end
    mu{iphase}(nc) = mu_last;        
    chk_mu(iphase) = mu{iphase}*p_out{pc_id(iphase)}(c_ind)'-g(iphase);  % Check that g = sum(mu*p)
end
% Display some results
disp(phases(pc_id))                                                      % Display stable phase assemblage
disp(phi')                                                               % Display volume fraction
disp(Cwt)                                                                % Display phase composition
disp(rhos)                                                               % Display solid density
disp(max(abs(chk_mu)))                                                   % Check consistency of chemical potential and g