clear,clf,addpath ../,addpath ../Solutions
% Fig6d
T = linspace(1000,1400,10) + 273.15;
P = linspace(0.01,1.2,11)*1e9;            
Cname = {'Si','Al'       ,'Fe'   ,'Mg', 'Ca', 'Na'   'Ti','O'};
Nsys =  [50.06  2*9.31    7.64   16.36 14.84  2*1.54  0.38];
noxy  = [2      3/2        1     1      1     1/2       2];
Nsys = [Nsys Nsys*noxy']; % someone measured rock composition of a pelite
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name    = {'Olivine','Clinopyroxene','Feldspar(I1)','Melt(H18)','Ilmenite'};
td        = init_thermo(phs_name,Cname,'solution_models_RE46');
[T2d,P2d] = ndgrid(T,P);
p = props_generate(td);
[g,Nphs,psc_id,p] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td); % compute Gibbs energy for all possible phases
LB  = zeros(1,size(g,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
parfor iTX = 1:length(T2d(:))
    alph{iTX} = linprog(g(:,iTX),[],[],Nphs,Nsys',LB); % The Gibbs energy minimization           
    disp(iTX/length(T2d(:)))
end
% Plotting
solv_tol = 60;
asm_id = zeros(length(T)*length(P),length(phs_name));
for iTX = 1:length(T2d(:))
    [alph_out,Npc,p_out,pc_id] = cluster_p(alph{iTX},Nphs,p,psc_id,solv_tol,phs_name);% Clustering algorithm
    asm_id(iTX,1:length(alph_out)) = pc_id; % Phase assemblage
end
figure(1),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name);
ylabel('P (GPa)'),xlabel('T(\circC)')