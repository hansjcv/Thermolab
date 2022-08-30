function [p_from_c_cons,icomp_indep,isite_indep] = comp2prop(comp_table,zt)
% Add order equations
% ord_eqn = zeros(length(ip_od),size(comp_table,2));
% for i_ord = 1:size(ord_eqn,1)
%     ord_eqn(i_ord,ip_od(i_ord)) = 1;
% end
p_constraints = [ones(1,size(comp_table,2)); comp_table; zt];
eq_id         = [1 2*ones(1,size(comp_table,1)) 3*ones(1,size(zt,1))];
% Find independent equations in p_constraints
eqns = p_constraints;
ieq_indep    = 1:size(p_constraints,1);
for icons = 1:size(p_constraints,1)
    itry       = 1:size(eqns,1);
    ieqn       = size(p_constraints,1) + 1 - icons;
    itry(ieqn) = [];
    if rank(eqns) == rank(eqns(itry,:))
        eqns           = eqns(itry,:);
        ieq_indep(end+1-icons) = -1;
    end
end
% ieq_indep(ieq_indep<0) = [];
% icomp_indep = ieq_indep(ieq_indep<=1+length(comp)&ieq_indep>1) - 1; % this assumes first equation is sum(p)=1
% isite_indep = ieq_indep(ieq_indep>length(comp)+1)-length(comp)-1; % this assumes first equation is sum(p)=1 and last equations are site fractions
icomp_indep = ieq_indep(eq_id==2&ieq_indep>0)-1;
isite_indep = ieq_indep(eq_id==3&ieq_indep>0)-sum(eq_id==2)-sum(eq_id==1);
ieq_indep(ieq_indep<0) = [];
p_from_c_cons = inv(p_constraints(ieq_indep,:));
end