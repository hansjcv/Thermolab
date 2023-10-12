clear,clf, addpath ../ ../EOS ../Solutions/ ../Utilities/
run_name = 'serp';
T      = linspace(200,900,10) + 273.15;
P      = linspace(0.1,4.0,11)*1e9;
solmod = 'solution_models_EF21';
ngrid    = 4; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 1;
Cname  = {'Si' ,'Fe'   ,  'Mg',   'H','O'  };
Nsys   = [34      10       38+1    62+2 147+2+0.01];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Antigorite','Brucite','Olivine','Orthopyroxene','Talc','Spinel','Lizardite','Fluid'};
td       = init_thermo(phs_name,Cname,solmod);
for i = 1:length(phs_name),td(i).dz(:) = 1/6;end
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