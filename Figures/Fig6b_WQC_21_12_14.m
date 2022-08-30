clear,clf,addpath ../
% Fig6b
T = linspace(700,1200,100) + 273.15;
P = 1.4e9;
Xsys  = linspace(0,1,101);            
X     = {'Si','Ca','C' ,'H','O'};
Nsys1 = [1     1    1e5   2    2+1+2e5+1];
Nsys2 = [1     1    1     2e5  2+1+2+1e5];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_id  = {'Fluid-CO2-H2O','wo,tc-ds55','q,tc-ds55','cc,tc-ds55'};
td            = init_thermo(phs_id,X,'solution_models_HP98');
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
xlabel('X_{CO_2}'),ylabel('T(\circC)')
hold on,
load data_Aranovich_2013_wo_cc
plot(data_14kb(:,1),data_14kb(:,2),'ok','MarkerSize',10)