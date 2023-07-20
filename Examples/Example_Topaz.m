clear,clf,addpath ../ ../Utilities/ ../Solutions/ ../EOS
% Example phase diagrams Fig 4
T = linspace(100,1050,20) + 273.15;
P = linspace(0.1,1,21)*1e9;  
Cname = {'Si','Al','K','F','H','O'};
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name  = {'q,tc-ds55','mic,tc-ds55','Topaz(DB)','Fluid(KF)','Crb,DB04','mu,tc-ds55','and,tc-ds55','HF(g),DB04'};
[T2d,P2d] = ndgrid(T,P); % (technical note: created 2D versions of T and P vector )
td = init_thermo(phs_name,Cname,'solution_models_HP98_Tpz');
td(3).dz(:) = 1/100;
p = props_generate(td);
[g,Nphs,pc_id] = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td,p); % compute Gibbs energy for all possible phases
Nsys = Nphs(:,pc_id==1) + Nphs(:,pc_id==2)  + 0.01*Nphs(:,pc_id==5) + [0 0 0 0 2 1]' + 0.5*Nphs(:,pc_id==7) + Nphs(:,pc_id==8);
LB  = zeros(1,size(g,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iPT = 1:length(T2d(:))
    alph{iPT} = linprog(g(:,iPT),[],[],Nphs,Nsys,LB); % The Gibbs energy minimization       
end
% Plotting
solv_tol = 20;
asm_id = zeros(length(T)*length(P),length(phs_name));
for iTX = 1:length(T2d(:))
    [alph_out,Npc,p_out,pc_id_out] = cluster_p(alph{iTX},Nphs,p,pc_id,solv_tol,phs_name);% Clustering algorithm
    asm_id(iTX,1:length(alph_out)) = pc_id_out; % Phase assemblage
end
figure(1),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name);
ylabel('P (GPa)'),xlabel('T(\circC)')