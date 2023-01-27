clear,clf
run_name = 'K32_2023_01_26_H18';
T      = linspace(400,900,10) + 273.15;
P      = linspace(0.1,2.0,11)*1e9;
solmod = 'solution_models_H18';
dz     = 1/4;
ngrid    = 1; % number of P-T grid refinements, each time, the P-T grid resolution is doubled.
eps_solv = 2;
Cname  = {'Si' ,'Al' , 'Cr',    'Ti'     ,'Fe'   ,'Mn',    'Mg',    'Ca',   'Na',    'K',    'H','O'  };
% wtOx   = [54.11  25.29  0        0.92   8.98    0*0.1     1.74    0.86     2.1    2.26     50];
% Oxname = {'SiO2','Al2O3','Cr2O3','TiO2',  'FeO'  ,'MnO', 'MgO',  'CaO',  'Na2O', 'K2O',  'H2O'};
wtOx   = [60.04  15.91  0        0.77     8.36     0.11   4.81    4.14    2.34    2.33     1.77];
% wtOx   = [40.73  2.52   0        0        7.65       0      36.07    2.7     0        0       50];
noxy   = [2      3      3        2        1       1         1        1       1        1       1       ];
ncat   = [1      2      2        1        1       1         1        1       2        2       2       ];
molmOx = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysOx = wtOx./molmOx;
Nsys   = NsysOx.*ncat;Nsys = [Nsys Nsys*(noxy./ncat)']; Nsys = Nsys/sum(Nsys);% someone measured rock composition of a pelite
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Quartz','Chlorite','Garnet','Spinel','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid',...
    'Antigorite','Brucite','Olivine','Orthopyroxene','Clinopyroxene','Amphibole','Talc','Magnesite','Epidote','Melt(H18)',...
    'Cordierite','Ilmenite','Rutile','Andalusite','Sillimanite','Kyanite','Fluid-H2O','Lime','Periclase','Corundum','Zoisite','Lawsonite'};
options.nref     = 1; % max number of iterations
    options.eps_dg   = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
    options.dz_tol   = 1e-14; % tolerance to stop iterations when z window becomes below this
    options.z_window = ones(size(phs_name))*0.085; % the window over which the refined grid is generated
    options.dz_fact  = ones(size(phs_name))*1.5; % the factor to determine new dz spacing, the larger, the more pseudocompounds
    options.ref_fact = 1.25; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
    options.disp_ref = 0; % show refinement graphically
    options.solver   = 0;
td       = init_thermo(phs_name,Cname,solmod);
c_exc = find(Nsys==0);
for i_sol = 1:length(phs_name)    
    exc_sol(i_sol) = sum(sum((td(i_sol).n_em(:,c_exc))>0,2)>0)==size(td(i_sol).n_em,1);    
end
Cname(Nsys==0) = [];%
Nsys(Nsys==0)  = [];
phs_name(exc_sol) = [];%
td       = init_thermo(phs_name,Cname,solmod);
for i_sol = 1:length(phs_name)
    td(i_sol).dz(:) = dz;
end
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