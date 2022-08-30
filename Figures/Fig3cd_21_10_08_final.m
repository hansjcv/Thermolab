clear,clf,addpath ../
% Input
T     = 1000+ 273.15;     % Temperature (K)
P     = 0.1e9;            % Pressure    (Pa)
phase = {'q,tc-ds55','per,tc-ds55','en,tc-ds55','fo,tc-ds55'};  % Phase name
g     = tl_gibbs_energy(T,P,phase);
Nsys = [1.2  % Si in system
        2.0  % Mg in system
        4.4];% O  in system
Nphs = [1 0 2 1  % Si in phases
        0 1 2 2  % Mg in phases
        2 1 6 4];% O  in phases
LB   = zeros(1,length(g));            % lower bounds for linprog minimization
alph = linprog(g,[],[],Nphs,Nsys,LB); % constrained Gibbs energy minimization
disp(phase(alph>0));
mol_fraction = alph/sum(alph);
Xphs = Nphs(1:2,:)./sum(Nphs(1:2,:));
gbar = g'./sum(Nphs(1:2,:));
gsys = g'*alph/sum(Nsys(1:2,:));
Xsys = Nsys(1:2,:)/sum(Nsys(1:2,:));
Xstb = Xphs(:,alph>0);
gstb = gbar(alph>0);
plot(Xphs(1,:),gbar,'ob',Xsys(1),gsys,'dr',Xstb(1,:),gstb,'*-g','LineWidth',1,'MarkerSize',10)
legend('Phases','System','Stable phases')
xlabel('SiO_2/(MgO+SiO_2) (mol)'),ylabel('g (J/mol)')
set(gca,'FontSize',10)
text(Xphs(1,:)+0.01,gbar+0.2e5,{'Qtz','Per','En','Fo'})
axis square