function mu = mu_fun(T,P,phases,td,p,g0,v0,rho_w,eps_w)
% Numerics
delc    = 1e-8;
g       = tl_gibbs_energy(T,P,phases,td,p,g0,v0,rho_w,eps_w);
mu      = cell(size(phases));
% % Chemical potentials
% for iphase = 1:length(phases)
%     Gphase  = g(iphase);
%     mu_last = Gphase;
%     c_ind   = find(p{iphase}~=0);                
%     [max_p,id_max_p]= max(p{iphase}(c_ind));
%     nc      = c_ind(id_max_p);
%     c_ind(id_max_p) = [];
%     if length(c_ind) >= 1
%         for ic = 1:length(c_ind)
%             p1                    = p;
%             p1{iphase}(c_ind(ic)) = p{iphase}(c_ind(ic)) + delc;
%             p1{iphase}(nc)        = p{iphase}(nc) - delc;
%             g_dc                  = tl_gibbs_energy(T,P,phases,td,p1,g0,v0,rho_w,eps_w);
%             Gphase1               = g_dc(iphase);
%             dGdc1                 = (Gphase1-Gphase)/delc;
%             mu_last               = mu_last - p{iphase}(c_ind(ic)).*dGdc1;
%             mu{iphase}(ic)        = dGdc1;
%         end
%         for ic = 1:length(c_ind)
%             mu{iphase}(ic) = mu{iphase}(ic) + mu_last;
%         end
%     end
%     mu{iphase}(nc) = mu_last;
% end
p_out = p;
pc_id = 1:length(phases);
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
%     chk_mu(iphase) = mu{iphase}*p_out{pc_id(iphase)}(c_ind)'-g(iphase)  % Check that g = sum(mu*p)
end