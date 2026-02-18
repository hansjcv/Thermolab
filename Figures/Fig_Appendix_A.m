clear
%% Amphibole solution model data
sitenames   = {'A',  'A' , 'A', 'M13','M13','M2', 'M2', 'M2','M2' , 'M2', 'M4','M4','M4','M4','T1','T1','V' ,'V'};
occupancy   = {'v',  'Na', 'K', 'Mg' ,'Fe' ,'Mg', 'Fe', 'Al','Fe3', 'Ti', 'Ca','Mg','Fe','Na','Si','Al','OH','O'};
st          = [
                1,    0,    0,   3,    0,    2,    0,    0,    0,    0,    2,   0,   0,   0,   4,   0,   2,   0;
                1,    0,    0,   3,    0,    0,    0,    2,    0,    0,    2,   0,   0,   0,   2,   2,   2,   0;
                0,    1,    0,   3,    0,    1,    0,    1,    0,    0,    2,   0,   0,   0,   2,   2,   2,   0;
                1,    0,    0,   3,    0,    0,    0,    2,    0,    0,    0,   0,   0,   2,   4,   0,   2,   0;
                1,    0,    0,   3,    0,    2,    0,    0,    0,    0,    0,   2,   0,   0,   4,   0,   2,   0;
                1,    0,    0,   0,    3,    0,    2,    0,    0,    0,    0,   0,   2,   0,   4,   0,   2,   0;
                1,    0,    0,   3,    0,    0,    2,    0,    0,    0,    0,   0,   2,   0,   4,   0,   2,   0;
                1,    0,    0,   0,    3,    2,    0,    0,    0,    0,    0,   0,   2,   0,   4,   0,   2,   0;
                1,    0,    0,   3,    0,    0,    0,    0,    2,    0,    0,   0,   0,   2,   4,   0,   2,   0;
                0,    0,    1,   3,    0,    1,    0,    1,    0,    0,    2,   0,   0,   0,   2,   2,   2,   0;
                1,    0,    0,   3,    0,    0,    0,    0,    0,    2,    2,   0,   0,   0,   2,   2,   0,   2];
Cname      = {'Si','Ti','Al','Ca','Fe','Mg','Na','K','H','O'};
comp_table = [
               8,   0,   0,   2,   0,   5,   0,   0,  2,  24;
               6,   0,   4,   2,   0,   3,   0,   0,  2,  24;
               6,   0,   3,   2,   0,   4,   1,   0,  2,  24;
               8,   0,   2,   0,   0,   3,   2,   0,  2,  24;
               8,   0,   0,   0,   0,   7,   0,   0,  2,  24;
               8,   0,   0,   0,   7,   0,   0,   0,  2,  24;
               8,   0,   0,   0,   4,   3,   0,   0,  2,  24;
               8,   0,   0,   0,   5,   2,   0,   0,  2,  24;
               8,   0,   0,   0,   2,   3,   2,   0,  2,  24;
               6,   0,   3,   2,   0,   4,   0,   1,  2,  24;
               6,   2,   2,   2,   0,   3,   0,   0,  0,  24]';
%% Find the site indices
sites   = unique(sitenames,'stable'); % List of all sites
site_id = zeros(1,length(sitenames'));
for i_site = 1:length(sites)
    site_id(strcmp(sitenames,sites(i_site))) = i_site; % index for the site (rather than a string with names)
end
%% Make site fraction table
max_st = zeros(size(st));
for i_site = 1:max(site_id)
    max_st(:,site_id==i_site) = repmat(sum(st(:,site_id==i_site),2),1,sum(site_id==i_site));% Find maximum site occupancy moles
end
zt = st./max_st;   zt(isnan(zt)) = 0;
%% Site fraction as function of proportion
z_equations          = [ones(1,size(zt,1)); zt'];% The proportions summing up to one, and z equations in one matrix
[z_eqns,ieq_z_indep] = indep_eqns(z_equations);  % Find independent site fraction equations
site_var             = ieq_z_indep(2:end)-1; % Find independent site variable index (this assumes first equation is sum(p)=1)
p_from_z_cons        = inv(z_eqns);              % Proportion from site fraction matrix
%% Bulk chemistry composition as function of proportion
p_equations          = [ones(1,size(comp_table,2)); comp_table; zt'];         % Augmented matrix of sum(p)=1, comp from p, and zt'
eq_id                = [1 2*ones(1,size(comp_table,1)) 3*ones(1,size(zt,2))]; % Identifier for equation type
[p_eqns,ieq_p_indep] = indep_eqns(p_equations);                               % Find independent equations
icomp_indep          = find(eq_id(ieq_p_indep)==2)-1;                         % Find independent compositional variables
isite_indep   = ieq_p_indep(eq_id(ieq_p_indep)==3)-sum(eq_id==2)-sum(eq_id==1); % Find site fraction variables determining order-disorder
p_from_c_cons = inv(p_eqns);                                                  % Proportion from composition matrix
%% Function for automated search for independent equations
function [eqns,ieq_indep] = indep_eqns(all_equations)
% Find independent equations in p_constraints
eqns = all_equations;
ieq_indep    = 1:size(all_equations,1);
for icons = 1:size(all_equations,1)
    itry       = 1:size(eqns,1);
    ieqn       = size(all_equations,1) + 1 - icons;
    itry(ieqn) = [];
    if rank(eqns) == rank(eqns(itry,:))
        eqns           = eqns(itry,:);
        ieq_indep(end+1-icons) = -1;
    end
end
ieq_indep(ieq_indep<0) = [];
end