function ieq_indep = indep_eqns(all_eqns)
eqns = all_eqns;
ieq_indep    = 1:size(all_eqns,1);
for icons = 1:size(all_eqns,1)
    itry       = 1:size(eqns,1);
    ieqn       = size(all_eqns,1) + 1 - icons;
    itry(ieqn) = [];
    if rank(eqns) == rank(eqns(itry,:))
        eqns           = eqns(itry,:);
        ieq_indep(end+1-icons) = -1;
    end
end
ieq_indep(ieq_indep<0) = [];