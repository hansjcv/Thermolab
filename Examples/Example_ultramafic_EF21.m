clear,clf, addpath ../ ../EOS ../Solutions/ ../Utilities/
run_name = 'EF21_red';
T      = linspace(300,850,20) + 273.15;
P      = linspace(0.5,4.0,21)*1e9;
solmod = 'solution_models_EF21';
ngrid    = 4; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 1;
Cname   = {'Si'   ,'Al'   ,'Mg', 'Fe',  'S', 'H','O'  };
A       = [1        0       0    0       0    0    2
           0        2       0    0       0    0    3
           0        0       1    0       0    0    1
           0        0       0    1       0    0    1
           0        0       0    2       0    0    3
           0        0       0    1       2    0    0
           0        0       0    0       0    2    1];
wtnames = {'SiO2' ,'Al2O3','MgO','FeO','Fe2O3','FeS2','H2O'};
% wtperc  = [34.17  1.81     33.44  2.23  4.33    0.12   23.90];
wtperc  = [34.31  1.81     33.58  5.88  0.30    0.12   24.00];
molm    = molmass_fun(wtnames,Cname,A);
Nsys_oxi = wtperc./molm';
Nsys      = Nsys_oxi*A;
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Antigorite','Brucite','Olivine','Orthopyroxene','Talc','Spinel','Chlorite','Garnet','Fluid','Pyrrhotite','Clinohumite','pyr,tc-ds633','anth,tc-ds633'};
td       = init_thermo(phs_name,Cname,solmod);
for i = 1:length(phs_name),td(i).dz(:) = 1/4;end
p         = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
for i_grid = 1:ngrid
    [T2d,P2d] = ndgrid(T,P);
    parfor iPT = 1:length(T2d(:))
        if refine_id(iPT) == 1
            [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td);
            disp(iPT/length(T2d(:)))
        end
    end
    % Postprocessing
    assemblage_id = zeros(length(T)*length(P),length(phs_name));
    for iPT = 1:length(T2d(:))
        [alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT}] = cluster_p(alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT},eps_solv,phs_name);
        assemblage_id(iPT,1:length(alph_all{iPT})) = pc_id_ref{iPT};
    end
    % PT grid refinement
    if i_grid < ngrid
        [ass_legend,ass_name,assemblage_id] = tl_psection(T,P,Cname,assemblage_id,phs_name);drawnow
        [alph_all, T,P,refine_id,Npc_all,pc_id_ref,p_ref] = refine_grid(T,P,alph_all,assemblage_id,p_ref,pc_id_ref,Npc_all);
    end
end
save(['linprog_run_' run_name],'-v7.3');
tl_psection(T-273.15,P/1e9,Cname,assemblage_id,phs_name,0,[0,0],8);