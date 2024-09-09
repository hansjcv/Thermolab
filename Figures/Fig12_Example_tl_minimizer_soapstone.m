clear,addpath ..\ ..\Solutions\ ..\Utilities\ ..\EOS
run_name = 'soapstone_2022_02_14_n500_nz15_full_0_5kb';
T = 300 + 273.15;
P = 0.5e8;
Cname    = {'Si'     ,'Al'      ,'Mg',     'Fe',    'Ca',     'C',      'H',      'O' ,'e'};
n_oxy    = [2          3/2        1         1         1         2         1/2             ];
Nsys0    = [0.7392    0.0451    1.0000    0.0644    0.0113    0.0232    1.3895   3.3630  0]; % Linna Atg + 2 perc Linna Do
% Cvar     = [linspace(1e-2,0.43989,500) linspace(0.4399,0.4488,500) linspace(0.44881,1.1,500)];
Cvar     = [linspace(1e-2,4.1,500)];
fluid    = 'Fluid-CO2-H2O';
phs_name = {fluid,'Chlorite','Antigorite','Clinopyroxene','Talc','Magnesite','Dolomite','Brucite','Olivine','per,tc-ds55','q,tc-ds55','ky,tc-ds55','sill,tc-ds55','and,tc-ds55'};
td        = init_thermo(phs_name,Cname,'solution_models_soapstone');
for iz = 1:length(td)
    td(iz).dz = ones(size(td(iz).dz))/15;
end
[T2d,P2d,X2d] = ndgrid(T,P,Cvar);
p         = props_generate(td);     % generate endmember proportions
options.nref     = 150; % max number of iterations
options.eps_dg   = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.dz_tol   = 1e-14; % tolerance to stop iterations when z window becomes below this
options.z_window = ones(size(phs_name))*0.2; % the window over which the refined grid is generated
options.dz_fact  = ones(size(phs_name))*3; % the factor to determine new dz spacing, the larger, the more pseudocompounds
options.ref_fact = 1.5; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
options.disp_ref = 0;  % display refinement iterations
options.solver = 0;
[rho_w,eps_di]   = water_props(T2d(:),P2d(:),phs_name,'ZD05','S14');       % Water properties
tic
% Minimization refinement
parfor iPT = 1:length(T2d(:))
    Nsys = Nsys0';
    Nsys(strcmp(Cname,'C')) = Nsys0(strcmp(Cname,'C')) + Cvar(iPT);
    Nsys(end-1) = n_oxy*Nsys(1:end-2);
    [alph_all{iPT},Npc_all{iPT},pc_id_all{iPT},p_all{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td,options,rho_w(iPT,:),eps_di(iPT,:));
    disp(iPT/length(T2d(:)))    
end
cpu = toc
save(['linprog_run_' run_name]);