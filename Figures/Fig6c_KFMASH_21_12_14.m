clear,clf
% Fig6c
T = linspace(400,950,20) + 273.15;
P = linspace(1,15,21)*1e8;            
Cname = {'Si','Al'       ,'Fe'   ,'Mg', 'K',    'H','O'};
noxy  = [2    3/2        1       1      1/2   1/2    ];
Nsys  = [68   2*12.488   10.529  4.761  2*3.819  2*50];Nsys = [Nsys Nsys*noxy']; % someone measured rock composition of a pelite
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
% phs_name = {'Chlorite','Garnet','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid','Cordierite','and,tc-ds55','sill,tc-ds55','ky,tc-ds55','H2O,tc-ds55','q,tc-ds55'};
% td         = init_thermo(phs_name,Cname,'solution_models_KFMASH');
phs_name = {'Melt(H18)','Chlorite','Garnet','Biotite','Muscovite','Staurolite','Feldspar(C1)','Chloritoid','Cordierite','and,tc-ds633','sill,tc-ds633','ky,tc-ds633','H2O,tc-ds633','q,tc-ds633'};
td         = init_thermo(phs_name,Cname,'solution_models_H18');
[T2d,P2d] = ndgrid(T,P); 
[g,Nphs,psc_id,p] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td); % compute Gibbs energy for all possible phases
LB  = zeros(1,size(g,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
parfor iTX = 1:length(T2d(:))
    alph{iTX} = linprog(g(:,iTX),[],[],Nphs,Nsys',LB); % The Gibbs energy minimization           
    disp(iTX/length(T2d(:)))
end
% Plotting
solv_tol = 20;
asm_id = zeros(length(T)*length(P),length(phs_name));
for iTX = 1:length(T2d(:))
    [alph_out,Npc,p_out,pc_id] = cluster_p(alph{iTX},Nphs,p,psc_id,solv_tol,phs_name);% Clustering algorithm
    asm_id(iTX,1:length(alph_out)) = pc_id; % Phase assemblage
end
figure(1),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name);
ylabel('P (GPa)'),xlabel('T(\circC)')