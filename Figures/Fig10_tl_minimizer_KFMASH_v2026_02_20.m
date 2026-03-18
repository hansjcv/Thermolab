clear,clf,addpath ../ ../Solutions ../Utilities
run_name = 'KFMASH_2026_02_20_80x81';
T     = linspace(400,700,80) + 273.15;
P     = linspace(0.1,1.5,81)*1e9;
Cname_oxy = {'SiO2','Al2O3'  ,'FeO'   ,'MgO'  , 'K2O' ,    'H2O'};
Cmol_oxy  = [68    , 12.488  , 10.529 , 4.761 , 3.819 ,       50];
[Nsys,Cname] = Oxidemol2Elementalmol(Cname_oxy,Cmol_oxy);
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Chlorite','Garnet','Biotite','Muscovite','Staurolite','Chloritoid',...
            'Cordierite','and,tc-ds55','sill,tc-ds55','ky,tc-ds55','H2O,tc-ds55','q,tc-ds55'};
td       = init_thermo(phs_name,Cname,'solution_models_KFMASH');
for i = 1:length(phs_name),td(i).nc(:) = 3;end
p        = props_generate(td);     % generate endmember proportions
[T2d,P2d] = ndgrid(T,P);
options.nref = 150;
% Minimization refinement
parfor iPT = 1:length(T2d(:))
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td,options);
    disp(iPT/length(T2d(:)))
end
save(['linprog_run_' run_name],'-v7.3');