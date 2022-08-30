clear,clf,addpath ../
% Fig6a
T = linspace(1200,2000,30) + 273.15;
P = 1e5;
Xsys  = linspace(0,1,31);            
X     = {'Si','Fe','Mg','O'};
Nsys1 = [1      0    2   4];
Nsys2 = [1      2    0   4];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_id        = {'Olivine','Melt'};
td            = init_thermo(phs_id,X,'solution_models_H18');
[X2d,T2d,P2d] = ndgrid(Xsys,T,P); % (technical note: created 2D versions of T and P vector )
[g,Nphs,psc_id,p] = tl_gibbs_energy(T2d(:),P2d(:),phs_id,td); % compute Gibbs energy for all possible phases
LB  = zeros(1,size(g,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iTX = 1:length(T2d(:))
    Nsys = Nsys1*X2d(iTX) + Nsys2*(1-X2d(iTX));
    alph{iTX} = linprog(g(:,iTX),[],[],Nphs,Nsys',LB); % The Gibbs energy minimization           
    disp(iTX/length(T2d(:)))
end
% Plotting
solv_tol = 20;
asm_id = zeros(length(Xsys)*length(T),length(phs_id));
for iTX = 1:length(T2d(:))
    [alph_out,Npc,p_out,pc_id] = cluster_p(alph{iTX},Nphs,p,psc_id,solv_tol,phs_id);% Clustering algorithm
    asm_id(iTX,1:length(alph_out)) = pc_id; % Phase assemblage
end
figure(1),
tl_psection(Xsys,T-273.15,X,asm_id,phs_id);
xlabel('X_{fo}'),ylabel('T(\circC)')
Tdata_liq = [1340 1450 1495];
Xdata_liq = [10   20   25]/100;
Tdata_sol = [1275 1317 1410 1440 1465 1495];
Xdata_sol = [20   30   46   52.5   55  59 ]/100;
hold on
plot(Xdata_liq,Tdata_liq,'ok','MarkerSize',10,'LineWidth',2)
plot(Xdata_sol,Tdata_sol,'ok','MarkerSize',10,'LineWidth',2)