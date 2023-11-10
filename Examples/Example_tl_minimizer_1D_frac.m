clear,clf,addpath ../ ../Utilities/ ../Solutions/ ../EOS
run_name = 'K32_2023_01_26_PTpath';
T        = linspace(600,850,20) + 273.15; % Temperature
dPdT     =  2e9/200;
P0       = 1e9;
P        = P0 + dPdT*(T-T(1));           % rho*g*z assumption
solmod = 'solution_models_HP98';
Oxname = {'SiO2','Al2O3','Cr2O3','TiO2',  'FeO'  ,'MnO', 'MgO',  'CaO',  'Na2O', 'K2O',  'H2O'};
wtOx   = [60.04  15.91  0        0.77     8.36     0.11   4.81    4.14    2.34    2.33     50];
dz     = 1/6;
solv_tol = 1;
Cname  = {'Si' ,'Al' , 'Cr',    'Ti'     ,'Fe'   ,'Mn',    'Mg',    'Ca',   'Na',    'K',    'H','O'  };
noxy   = [2      3      3        2        1       1         1        1       1        1       1       ];
ncat   = [1      2      2        1        1       1         1        1       2        2       2       ];
molmOx = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysOx = wtOx./molmOx;
Nsys   = NsysOx.*ncat;Nsys = [Nsys Nsys*(noxy./ncat)']; Nsys = Nsys/sum(Nsys);% someone measured rock composition of a pelite
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Chlorite','Garnet','Spinel','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid',...
    'Antigorite','Brucite','Amphibole','Olivine','Orthopyroxene','Clinopyroxene','Talc','Magnesite','Epidote',...
    'Cordierite','Ilmenite','ru,tc-ds55','and,tc-ds55','sill,tc-ds55','ky,tc-ds55','H2O,tc-ds55','lime,tc-ds55','per,tc-ds55','cor,tc-ds55','zo,tc-ds55','q,tc-ds55'}; 
td       = init_thermo(phs_name,Cname,solmod);
c_exc = find(Nsys==0);
for i_sol = 1:length(phs_name)    
    exc_sol(i_sol) = sum(sum((td(i_sol).n_em(:,c_exc))>0,2)>0)==size(td(i_sol).n_em,1);    
end
Cname(Nsys==0)    = [];%
Nsys(Nsys==0)     = [];
phs_name(exc_sol) = [];%
td       = init_thermo(phs_name,Cname,solmod);
for i_sol = 1:length(phs_name)
    td(i_sol).dz(:) = dz;
end
p        = props_generate(td);     % generate pseudocompounds
tic
Nsyst = Nsys;
% Minimization refinement
for iPT = 1:length(T(:))
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T(iPT),P(iPT),Nsys,phs_name,p,td);
    [alph,Npc,p_out,pc_id] = cluster_p(alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT},solv_tol,phs_name);    % Compute exsolved, or numerically equivalent phases 
    if sum(pc_id==3)>0
        Nsys =  Nsys - Npc(:,pc_id==3)'*alph(pc_id==3);     
    end
    disp(iPT/length(T(:)))
end
cpu = toc
save(['linprog_run_' run_name],'-v7.3');