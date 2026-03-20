clear,clf,addpath ../ ../Solutions ../Utilities ../EOS
T = 1700 + 273.15;%800 + 273.15;
P = 40e8;%26.7424*1e8;         
Cname_oxy = {'SiO2','Al2O3'  ,'CaO'   ,'MgO'  , 'FeO' ,'K2O','Na2O', 'TiO2','Cr2O3'};
nO        =  0.096;
Cmol_oxy  = [38.49 1.78 2.82 50.57 5.89 0.010 0.250 0.100 0.096 0.109];
[Nsys,Cname] = Oxidemol2Elementalmol(Cname_oxy,Cmol_oxy);
Nsys(strcmp(Cname,'O')) = Nsys(strcmp(Cname,'O')) + nO; %
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name    = {'Olivine','Orthopyroxene','Clinopyroxene','Feldspar','Spinel','Garnet','Melt(G25)','ru,tc-ds636','q,tc-ds636'};
td        = init_thermo(phs_name,Cname,'Igneous');
for i_p = 1:length(phs_name),td(i_p).nc(:) = 3;end
p        = props_generate(td);     % generate endmember proportions
options.disp_ref    = 1;
options.nref        = 150;
options.eps_dg      = 1e-3;
options.fsolve      = 1;
options.use_pgrid   = 1;
options.TPDmin      = 1;
% Minimization refinement
[alph_all,Npc_all,pc_id_ref,p_ref,g_min,alph_fsolve,Npc_fsolve,pc_id_fsolve,p_fsolve] = tl_minimizer(T,P,Nsys,phs_name,p,td,options);
molm = molmass_fun(Cname);
iPT = 1;
solv_tol = 1;
print_reactions = 1;
if ~isempty(alph_fsolve) && sum(alph_fsolve<0)==0 && isreal(alph_fsolve)
    Stablephases    = phs_name(pc_id_fsolve)';
    td_stable       = td(pc_id_fsolve);
    pc_id_fsolve    = 1:numel(pc_id_fsolve);
    [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T,P,td_stable,alph_fsolve,Npc_fsolve,molm,p_fsolve,pc_id_fsolve,Stablephases,solv_tol,'PS94','S14');
    postprocess_reactions(T,P,td_stable,pc_id,p_out,print_reactions);
else
    Stablephases    = phs_name(pc_id)';
    [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T,P,td,alph_all,Npc_all,molm,p_ref,pc_id_ref,phs_name,solv_tol,'PS94','S14');
    postprocess_reactions(T,P,td,pc_id,p_out,print_reactions);
end
%% Here some more human readable meaning of the output from the postprocess function
Volumefractions = phi; 
Weightfractions = phiw;
Molfractions    = molphase/sum(molphase);
Densities       = rho; % in kg/m3
GibbsEnergies   = g;  % in J/mol 
% Make table
EquilibriumResult = [table(Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies,'RowNames',Stablephases) array2table(sfu','VariableNames',Cname)];
% Print in command window
EquilibriumResult