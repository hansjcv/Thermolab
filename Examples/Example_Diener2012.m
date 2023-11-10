clear,clf, addpath ../ ../EOS ../Solutions/ ../Utilities/
run_name = 'mafic2';
T      = linspace(550,675,14) + 273.15;
P      = linspace(0.2,0.6,15)*1e9;
solmod = 'Metabasite_HP98';
ngrid    = 1; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 2;
Cname  = {'Si' ,'Ti',  'Al'     , 'Fe', 'Mg', 'Ca', 'Na',  'K',  'H','O'  };
Nsys   = [64.87  0.33   5.18*2  18.64   6.96  2.4    0.29*2 0.08*2  50   1.24];
% Nsys   = [34     0     0          0     48    0    0       0       62   0];
noxy   = [2       2     3        1     1     1     1        1    1 ];
ncat   = [1        1    2        1      1    1     2        2  2 ];
Nsys = [Nsys(1:end-1) Nsys(1:end-1)*(noxy./ncat)'+Nsys(end)]; Nsys = Nsys/sum(Nsys);
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Clino-amphibole','Ortho-amphibole','Clinopyroxene','Chlorite','Cordierite','Biotite','Muscovite','Garnet','Feldspar(C1)','Talc','Spinel','Ilmenite','H2O,tc-ds55','per,tc-ds55','q,tc-ds55'};
td       = init_thermo(phs_name,Cname,solmod);
options.nref     = 150; % max number of iterations
options.eps_dg   = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.dz_tol   = 1e-14; % tolerance to stop iterations when z window becomes below this
options.z_window = ones(size(phs_name))*0.085; % the window over which the refined grid is generated
options.dz_fact  = ones(size(phs_name))*1.5; % the factor to determine new dz spacing, the larger, the more pseudocompounds
options.ref_fact = 1.25; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
options.disp_ref = 0; % show refinement graphically
options.solver   = 0;
for i = 1:length(phs_name),td(i).dz(:) = 1/3;end
p         = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
for i_grid = 1:ngrid
    [T2d,P2d] = ndgrid(T,P);
    parfor iPT = 1:length(T2d(:))
        if refine_id(iPT) == 1
            [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td,options);
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