clear,clf
run_name = 'KFMASH_2022_10_22_ngrid2';
T     = linspace(400,750,6) + 273.15;
P     = linspace(0.1,1.5,5)*1e9;
Cname = {'Si','Al'       ,'Fe'   ,'Mg', 'K',    'H','O'};
noxy  = [2    3/2        1       1      1/2   1/2    ];
Nsys  = [68   2*12.488   10.529  4.761  2*3.819  2*50];Nsys = [Nsys Nsys*noxy']; % someone measured rock composition of a pelite
molm  = molmass_fun(Cname); % get molar mass of the components
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Chlorite','Garnet','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid',...
    'Cordierite','and,tc-ds55','sill,tc-ds55','ky,tc-ds55','H2O,tc-ds55','q,tc-ds55'};
td       = init_thermo(phs_name,Cname,'solution_models_HP98');
for i_sol = 1:length(phs_name)
    td(i_sol).dz(:) = 1/6;
end
p        = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
ngrid    = 4; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 100;
for i_grid = 1:ngrid
    [T2d,P2d] = ndgrid(T,P);
    for iPT = 1:length(T2d(:))
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