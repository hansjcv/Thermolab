clear,clf,addpath ../ ../Utilities/ ../Solutions/ ../EOS
run_name = 'example_PTpoint';
T      = 800 + 273.15; %K
P      = 4e9; % Pa
solmod = 'solution_models_H18'; % solution model file name. See Solutions subfolder for more
Cname  = {'Si' ,'Al', 'Fe'  'Mg',    'H','O'  }; % System components
Nsys   = [34     2     2    46        62  150.001]; % system composition in elemental mole.
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Quartz','Chlorite','Garnet','Spinel','Antigorite','Brucite','Olivine','Orthopyroxene','Amphibole','Talc','per,tc-ds633','H2O,tc-ds633'};
td       = init_thermo(phs_name,Cname,solmod); % initialize thermodynamic data
p        = props_generate(td);                 % generate pseudocompounds
% One minimization
[alph,Npc,pc_id,p,g_min] = tl_minimizer(T,P,Nsys,phs_name,p,td);
% Postprocess
molm = molmass_fun(Cname); % get molarmasses of the Cname;
solv_tol = 1; % tolerance for recognizing solvi
[pc_id,phi,Cwt,sfu,rho,mu,prop,phiw,g,molphase] = postprocess_fun(T,P,td,alph,Npc,molm,p,pc_id,phs_name,solv_tol); % postprocess minimization result
% prepare table variables for display only. 
% Here some more human readable meaning of the output from the postprocess function
Stablephases    = phs_name(pc_id)';
Volumefractions = phi; 
Weightfractions = phiw;
Molfractions    = molphase/sum(molphase);
Densities       = rho; % in kg/m3
GibbsEnergies   = g;  % in J/mol 
% Make table
EquilibriumResult = [table(Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies,'RowNames',Stablephases) array2table(sfu','VariableNames',Cname)];
% Print in command window
EquilibriumResult