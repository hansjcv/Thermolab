clear,clf, addpath ../ ../Solutions/ ../EOS
T = 400 + 273.15;
P = 1e9;
dP = 1e5;
phs_name = {'atg,tc-ds633','br,tc-ds633','fo,tc-ds633','H2O,tc-ds633'};
Cname = {'Si','Mg','H','O'};
td = init_thermo(phs_name,Cname);
molm       = molmass_fun(Cname);
g          = tl_gibbs_energy(T,P   ,phs_name);
g_plus_dP  = tl_gibbs_energy(T,P+dP,phs_name);
g_minus_dP = tl_gibbs_energy(T,P-dP,phs_name);
molweight  = cell2mat({td.n_em}')*molm;
V    = (g_plus_dP - g)/dP;
rho  = molweight./V;
V2   = (g - g_minus_dP)/dP;
rho2 = molweight./V2;
beta = -1./V.*(V-V2)/(dP);