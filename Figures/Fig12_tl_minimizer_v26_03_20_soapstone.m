clear,addpath ..\ ..\Solutions\ ..\Utilities\ ..\EOS
run_name = 'soapstone_2026_03_20_full_0_5kb';
T = 300 + 273.15;
P = 0.5e8;
Cname    = {'Si'     ,'Al'      ,'Mg',     'Fe',    'Ca',     'C',      'H',      'O' };
n_oxy    = [2          3/2        1         1         1         2         1/2             ];
Nsys0    = [0.7392    0.0451    1.0000    0.0644    0.0113    0.0232    1.3895   3.3630]; % Linna Atg + 2 perc Linna Do
Cvar     = [linspace(1e-2,0.463,500) linspace(0.4631,0.836,500) linspace(0.8361,4.1,500)];
nO       = 1e-3;
fluid    = 'Fluid-CO2-H2O';
phs_name = {fluid,'Chlorite','Antigorite','Clinopyroxene','Talc','Magnesite','Dolomite','Brucite','Olivine','mt,tc-ds55','per,tc-ds55','q,tc-ds55','ky,tc-ds55','sill,tc-ds55','and,tc-ds55'};
td        = init_thermo(phs_name,Cname,'solution_models_soapstone');
[T2d,P2d,X2d] = ndgrid(T,P,Cvar);
p         = props_generate(td);     % generate endmember proportions
options.nref      = 150; % max number of iterations
options.eps_dg    = 1e-3; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.TPDmin    = 1; % use tangent plane distance
options.fsolve    = 1;
[rho_w,eps_di]    = water_props(T2d(:),P2d(:),phs_name,'CORK','S14');       % Water properties
tic
% Minimization refinement
parfor iPT = 1:length(T2d(:))
    Nsys = Nsys0';
    Nsys(strcmp(Cname,'C')) = Nsys0(strcmp(Cname,'C')) + Cvar(iPT);
    Nsys(end) = n_oxy*Nsys(1:end-1)+nO;
    [alph_all{iPT},Npc_all{iPT},pc_id_all{iPT},p_all{iPT},gmin_all{iPT},alph_fsolve{iPT},Npc_fsolve{iPT},pc_id_fsolve{iPT},p_fsolve{iPT},chk_dg_ref(iPT),chk_Nsys_ref(iPT)] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys',phs_name,p,td,options,rho_w(iPT,:),eps_di(iPT,:));
    disp(iPT/numel(T2d))    
end
cpu = toc
save(['linprog_run_' run_name]);