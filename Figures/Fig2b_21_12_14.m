clear,clf,addpath ../
T  = (100:10:1000) + 273.15;
P  = [5, 10, 15, 20]*1e8;
linesym  = {'-b','-r','-g','-c'};
pointsym = {'ob','or','og','oc'};
legnames = {};
R  = 8.3144;
phases = {'aqSi,tc-ds55','q,tc-ds55'};
[T2d,P2d] = ndgrid(T,P);
td = init_thermo(phases);
p  = props_generate(td);
rho_w     = rho_H2O(T,P,'CORK');
eps_di    = eps_H2O(T,P,rho_w,'JN91');
[g0,v0]   = tl_g0(T2d(:),P2d(:),td,rho_w(:),eps_di(:));
[G,Nphs] = tl_gibbs_energy(T2d(:),P2d(:),phases,td,p,g0,v0,rho_w(:),eps_di(:));
v  = null(Nphs);v = v/v(1); % Find stoichiometric coefficients of reaction
dG = G'*v; 
dG = reshape(dG,length(T),length(P));
for iP = 1:length(P)
    % The reaction approach
    K       = exp(-dG(:,iP)/R./T');
    % Plotting    
    legnames = [legnames ['P = ' num2str(P(iP)/1e9) ' GPa']];
    h(iP) = plot(T-273.15,log10(K),linesym{iP},'LineWidth',1);
    xlabel('T(\circC)'),ylabel('log(mSiO_2)')
    % legend('q = coe')
    [data,txt] = xlsread('Manning94_data');
    T_data = data(:,1);
    P_data = data(:,2);
    mSiO2_data = data(:,3);
    wSiO2_data = data(:,4);
    hold on
    plot(T_data(P_data == P(iP)/1e8),log10(mSiO2_data(P_data == P(iP)/1e8)),pointsym{iP},'LineWidth',1)
end
legend(h,legnames,'Location','NorthWest')
set(gca,'FontSize',10)
axis square,axis tight
% print('Fig2b','-depsc','-r600');