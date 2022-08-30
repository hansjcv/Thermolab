clear,clf,addpath ../
% Example phase diagrams Fig 4
T = linspace(400,750,20) + 273.15;
P = linspace(0.1,1,21)*1e9;  
X = {'Si','Al','O'};
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phases  = {'and,tc-ds55','sill,tc-ds55','ky,tc-ds55'};
[T2d,P2d] = ndgrid(T,P); % (technical note: created 2D versions of T and P vector )
[g,Nphs] = tl_gibbs_energy(T2d(:),P2d(:),phases); % compute Gibbs energy for all possible phases
Nsys = Nphs(:,1);
LB  = zeros(1,size(g,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iPT = 1:length(T2d(:))
    alph(iPT,:) = linprog(g(:,iPT),[],[],Nphs,Nsys,LB); % The Gibbs energy minimization       
end
% Plotting
asm_id = zeros(length(T)*length(P),length(phases));
for i = 1:size(asm_id,1)
    asm_id(i,1:length(find(alph(i,:)>0))) = find(alph(i,:)>0);
end
figure(1),
tl_psection(T-273.15,P/1e9,X,asm_id,phases);
xlabel('T (\circC)'),ylabel('P(GPa)')
set(gca,'FontSize',12)