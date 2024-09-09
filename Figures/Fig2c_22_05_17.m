clear,clf,addpath ../ ../EOS
T   = [400 500 600]+273.15;
P   = linspace(0.001,3.5,500)*1e9;
linesym  = {'-b','-r','-g','-c'};
pointsym = {'ob','or','og','oc'};
legnames = {};
R  = 8.3144;
[T2d,P2d] = ndgrid(T,P);
for ipl = 1:2
    if ipl == 1
        phases = {'KCl0,Miron','K+,Miron','Cl-,Miron'};
        sym = '-';
    elseif ipl == 2
        phases = {'KCl,aq,supcrt','K+,supcrt','Cl-,supcrt'};
        sym = '--';
    end
    td        = init_thermo(phases);
    p         = props_generate(td);
    [rho_w,eps_di] = water_props(T2d(:),P2d(:),phases,'JN91','JN91');    
    [g0,v0]        = tl_g0(T2d(:),P2d(:),td,rho_w,eps_di);
    [G,Nphs]       = tl_gibbs_energy(T2d(:),P2d(:),phases,td,p,g0,v0,rho_w,eps_di);
    v  = null(Nphs);v = v/v(1); % Find stoichiometric coefficients of reaction
    dG = G'*v;
    lnK   = -dG./R./T2d(:);
    K     = exp(lnK);
    logK  = log10(K);
    plot(P/1e8',reshape(logK,length(T),length(P))',sym),
    set(gca,'ColorOrder',eye(3));
    axis([0 3.75 0 4.5])
    title('K^+ + Cl^- = KCl^0')
    xlabel('P(kb)'),ylabel('logK')
    % Plot Experimental data points
    load data_Miron16_KCl
    data_400C = [x_data(1:11) y_data(1:11)];
    data_500C = [x_data(12:21) y_data(12:21)];
    data_600C = [x_data(22:end) y_data(22:end)];
    % Plotting
    hold on
end
plot(data_400C(:,1),data_400C(:,2),'ro',data_500C(:,1),data_500C(:,2),'gd',data_600C(:,1),data_600C(:,2),'b*')
legend([num2str([T'-273.15])])