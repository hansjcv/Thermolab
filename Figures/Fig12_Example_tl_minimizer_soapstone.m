clear,addpath ..\ ..\Solutions\ ..\Utilities\ ..\EOS
run_name = 'soapstone_2022_02_14_n150_nz15_full_0_5kb';
T = 300 + 273.15;
P = 0.5e8;
Cname    = {'Si'     ,'Al'      ,'Mg',     'Fe',    'Ca',     'C',      'H',      'O' ,'e'};
n_oxy    = [2          3/2        1         1         1         2         1/2             ];
Nsys0    = [0.7392    0.0451    1.0000    0.0644    0.0113    0.0232    1.3895   3.3630  0]; % Linna Atg + 2 perc Linna Do
% Cvar     = [linspace(1e-2,0.43989,500) linspace(0.4399,0.4488,500) linspace(0.44881,1.1,500)];
Cvar     = [linspace(1e-2,4.1,150)];
fluid    = 'Fluid-CO2-H2O';
phs_name = {fluid,'Chlorite','Antigorite','Clinopyroxene','Talc','Magnesite','Dolomite','Brucite','Olivine','per,tc-ds55','q,tc-ds55','ky,tc-ds55','sill,tc-ds55','and,tc-ds55'};
td        = init_thermo(phs_name,Cname,'solution_models_soapstone');
for iz = 1:length(td)
    td(iz).nc = 3;
end
[T2d,P2d,X2d] = ndgrid(T,P,Cvar);
p         = props_generate(td);     % generate endmember proportions
options.nref      = 150; % max number of iterations
options.eps_dg    = 1e-5; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.disp_ref  = 0;  % display refinement iterations
options.disp_npc  = 0;  % display refinement pseudocompounds
options.TPDmin    = 1; % use tangent plane distance
options.use_pgrid = 1; % use proportions pseudocompound grid
[rho_w,eps_di]    = water_props(T2d(:),P2d(:),phs_name,'ZD05','S14');       % Water properties
tic
% Minimization refinement
parfor iPT = 1:length(T2d(:))
    Nsys = Nsys0';
    Nsys(strcmp(Cname,'C')) = Nsys0(strcmp(Cname,'C')) + Cvar(iPT);
    Nsys(end-1) = n_oxy*Nsys(1:end-2);
    [alph_all{iPT},Npc_all{iPT},pc_id_all{iPT},p_all{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td,options,rho_w(iPT,:),eps_di(iPT,:));
    disp(iPT/numel(T2d))    
end
cpu = toc
save(['linprog_run_' run_name]);