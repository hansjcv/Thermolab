clear,clf,addpath ../ ../Solutions ../Utilities ../EOS
T = 1000+273.15;%1252.476+273.15;%1080 + 273.15;
P = 1*1e8;            
Cname_oxy = {'SiO2','Al2O3'  ,'CaO'   ,'MgO'  , 'FeO' ,'K2O','Na2O', 'TiO2','Cr2O3'};
nO        =  0.35;
Cmol_oxy  = [50.72 9.16 15.21 16.25 7.06 0.01 1.47 0.39  0.01];
[Nsys,Cname] = Oxidemol2Elementalmol(Cname_oxy,Cmol_oxy);
Nsys(strcmp(Cname,'O')) = Nsys(strcmp(Cname,'O')) + nO; %
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name    = {'Olivine','Orthopyroxene','Clinopyroxene','Feldspar','Ilmenite','Melt(G25)','Spinel'};
td        = init_thermo(phs_name,Cname,'Igneous');
for i_p = 1:length(phs_name),td(i_p).nc(:) = 3;end
p        = props_generate(td);     % generate endmember proportions
td(7).nc = 6;
options.disp_ref  = 1;
options.disp_npc  = 0;
options.nref      = 150;
options.eps_dg    = 1e-8;
options.TPDmin    = 1;
options.gridmin   = 1;
options.fsolve    = 1;
options.solv_tol  = 1;
%% Minimization refinement
[alph_all,Npc_all,pc_id_ref,p_ref,g_min,alph_fsolve,Npc_fsolve,pc_id_fsolve,p_fsolve] = tl_minimizer(T,P,Nsys,phs_name,p,td,options);
%% Postprocess
molm = molmass_fun(Cname);
solv_tol = 1;
print_reactions = 1;
if ~isempty(alph_fsolve) && sum(alph_fsolve<0)==0 && isreal(alph_fsolve)
    Stablephases    = phs_name(pc_id_fsolve)';
    td_stable       = td(pc_id_fsolve);
    pc_id_fsolve    = 1:numel(pc_id_fsolve);
    [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T,P,td_stable,alph_fsolve,Npc_fsolve,molm,p_fsolve,pc_id_fsolve,Stablephases,solv_tol,'PS94','S14');
    postprocess_reactions(T,P,td_stable,pc_id,p_out,print_reactions);
else    
    [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T,P,td,alph_all,Npc_all,molm,p_ref,pc_id_ref,phs_name,solv_tol,'PS94','S14');
    postprocess_reactions(T,P,td,pc_id,p_out,print_reactions);
    Stablephases    = phs_name(pc_id)';
end
% Here some more human readable meaning of the output from the postprocess function
Volumefractions = phi; 
Weightfractions = phiw;
Molfractions    = molphase/sum(molphase);
Densities       = rho; % in kg/m3
GibbsEnergies   = g;  % in J/mol 
% Make table
EquilibriumResult = [table(Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies,'RowNames',Stablephases) array2table(sfu','VariableNames',Cname)];
% Print in command window
EquilibriumResult