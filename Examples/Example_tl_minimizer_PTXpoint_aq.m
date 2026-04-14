clear,clf,addpath ../ ../Utilities/ ../Solutions/ ../EOS
run_name = 'example_PTpoint';
T      = 600 + 273.15; %K
P      = 2e9; % Pa
solmod = 'Ultramafic'; % solution model file name. See Solutions subfolder for more
Cname  = {'Si' ,'Al', 'Fe'  'Mg',    'H','O'  'e'}; % System components
Nsys   = [34     2     2    46+1        62+2  150+0.001+2 0]; % system composition in elemental mole.
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
% phs_name = {'Chlorite','Garnet','Spinel','Antigorite','Brucite','Olivine','Orthopyroxene','anth,tc-ds633','Talc','q,tc-ds633','per,tc-ds633','cor,tc-ds633','H2O,tc-ds633'};
phs_name = {'Aq-fluid','Chlorite','Garnet','Spinel','Antigorite','Brucite','Olivine','Orthopyroxene','anth,tc-ds633','Talc','q,tc-ds633','per,tc-ds633'};
spc_name = {'SiO2,aq,DEW','Si2O4,aq,DEW','MgO,aq,DEW','MgOH+,DEW','Mg+2,DEW','Fe(HSiO3)+,DEW','AlO2(SiO2)-,DEW','AlO2-,DEW','HAlO2,aq,DEW','H+,DEW','OH-,DEW','H2O,tc-ds633'};
td       = init_thermo(phs_name,Cname,solmod,1,spc_name,1,0.5); % initialize thermodynamic data
td(1).subdtype(:) = 4;
for ip = 1:length(phs_name),td(ip).nc(:) = 3;end
p        = props_generate(td);                 % generate pseudocompounds
options.disp_ref   = 1;
options.TPDmin     = 1;
options.eps_dg     = 1e-4;
options.fsolve     = 1;
options.use_pgrid  = 1;
options.x0_old     = 1;
options.show_react = 1;
% One minimization
[alph,Npc,pc_id,p,g_min,alph_fsolve,Npc_fsolve,pc_id_fsolve,p_fsolve,dg_fsolve,Nsys_chk] = tl_minimizer(T,P,Nsys,phs_name,p,td,options);
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
[B,I] = sort((prop{1}'),'descend');
AqSpeciesList = td(1).p_name(I);
SpeciesMolefraction =  B;
SpeciesTable = table(SpeciesMolefraction,'RowNames',AqSpeciesList);
% Make table
EquilibriumResult = [table(Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies,'RowNames',Stablephases) array2table(sfu','VariableNames',Cname)];
% Print in command window
SpeciesTable
EquilibriumResult
dg_fsolve
Nsys_chk