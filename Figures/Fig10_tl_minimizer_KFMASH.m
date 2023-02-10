clear,clf,addpath ../,addpath ../Solutions
run_name = 'KFMASH_2022_03_06_10x11';
T     = linspace(400,700,10) + 273.15;
P     = linspace(0.1,1.5,11)*1e9;
Cname = {'Si','Al'       ,'Fe'   ,'Mg', 'K',    'H','O'};
noxy  = [2    3/2        1       1      1/2   1/2    ];
Nsys  = [68   2*12.488   10.529  4.761  2*3.819  2*50];Nsys = [Nsys Nsys*noxy']; % someone measured rock composition of a pelite
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Chlorite','Garnet','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid',...
            'Cordierite','and,tc-ds55','sill,tc-ds55','ky,tc-ds55','H2O,tc-ds55','q,tc-ds55'};
td       = init_thermo(phs_name,Cname,'solution_models_KFMASH');
% phs_name = {'Chlorite','Garnet','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid',...
%             'Cordierite','and,tc-ds633','sill,tc-ds633','ky,tc-ds633','H2O,tc-ds633','q,tc-ds633'};
% td       = init_thermo(phs_name,Cname,'solution_models_H18_simpler');
p        = props_generate(td);     % generate endmember proportions
[T2d,P2d] = ndgrid(T,P);
% Minimization refinement
parfor iPT = 1:length(T2d(:))
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td);
    disp(iPT/length(T2d(:)))
end
save(['linprog_run_' run_name],'-v7.3');