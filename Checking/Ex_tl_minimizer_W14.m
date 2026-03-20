clear,clf,addpath ../ ../Solutions ../Utilities ../EOS
T            = 750+273.15;%T2d(err_PT(end))%;550 + 273.15;
P            = 5e8;%P2d(err_PT(end))%1e8;
Cname_Oxide  = {'SiO2','Al2O3','CaO','MgO','FeO','K2O','Na2O','TiO2','H2O'};
molOx        = [69.198 8.836 0.276 3.594 7.039 2.833 0.563 0.616 6.516];Osys         = 0.528;
[Nsys,Cname] = Oxidemol2Elementalmol(Cname_Oxide,molOx); Nsys(strcmp(Cname,'O')) = Nsys(strcmp(Cname,'O')) + Osys;Nsys = Nsys/sum(Nsys);
phs_name    = {'Feldspar','Ilmenite','Garnet','Staurolite','Epidote','Muscovite','Margarite','Paragonite','Biotite','Chlorite','Chloritoid','Cordierite','Magnetite','Orthopyroxene','Melt(W14)','ru,tc-ds62','sph,tc-ds62','cor,tc-ds62','and,tc-ds62','ky,tc-ds62','sill,tc-ds62','q,tc-ds62','H2O,tc-ds62'};
td           = init_thermo(phs_name,Cname,'Metapelite');
td(strcmp(phs_name,'Melt(W14)')).subdtype = [4 4 4 4 0 4 4 0 0 0];
for i_p = 1:length(phs_name),td(i_p).nc(:) = 6;end
p        = props_generate(td);     % generate endmember proportions
options.disp_ref  = 1;
options.disp_npc  = 0;
options.nref      = 150;
options.eps_dg    = 1e-3;
options.TPDmin    = 1;
options.fsolve    = 1;
options.use_pgrid = 1;
%% Minimization refinement
[alph_all,Npc_all,pc_id_ref,p_ref,g_min,alph_fsolve,Npc_fsolve,pc_id_fsolve,p_fsolve] = tl_minimizer(T,P,Nsys,phs_name,p,td,options);
%% Postprocess
molm = molmass_fun(Cname);
solv_tol = 0.1;
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
EquilibriumResult = [table(Stablephases,Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies) array2table(sfu','VariableNames',Cname)];
% Print in command window
EquilibriumResult