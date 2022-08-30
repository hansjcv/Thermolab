function [p_from_z_cons,site_var,cons] = sites2prop(zt)
z_constraints = [ones(1,size(zt,1)); zt'];
% Find independent equations in p_constraints
eqns = z_constraints;
ieq_indep    = 1:size(z_constraints,1);
for icons = 1:size(z_constraints,1)
    itry       = 1:size(eqns,1);
    ieqn       = size(z_constraints,1) + 1 - icons;
    itry(ieqn) = [];
    if rank(eqns) == rank(eqns(itry,:))
        eqns           = eqns(itry,:);
        ieq_indep(end+1-icons) = -1;
    end
end
ieq_indep(ieq_indep<0) = [];
site_var             = ieq_indep(2:end) - 1; % Find independent site variable index (this assumes first equation is sum(p)=1)
p_from_z_cons        = inv(eqns);              % Proportion from site fraction matrix
if length(ieq_indep) == 1
    cons = [];
else
    cons = 1;
end
end