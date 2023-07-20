clear,addpath ../ ../Utilities/ ../Solutions/ ../Figures/ ../EOS
T     = linspace(400,900,20) + 273.15;
P     = 1e9;
phase = {'q,tc-ds633','H2O,tc-ds633','H+,DEW','OH-,DEW','SiO2,aq,DEW','HSiO3-,DEW','Si2O4,aq,DEW','Si3O6,aq,DEW','Si2O4OH-,DEW'};
m_id  = [0             1            2           2            2           2          2              2             2];
chrg  = [0             0            1          -1            0          -1          0              0            -1];
Cname = {'Si','H','O','e'};
td    = init_thermo(phase,Cname);
[T2d,P2d] = ndgrid(T,P);
[g0,Nphs] = tl_gibbs_energy(T2d(:),P2d(:),phase,td);
[rho_w,eps_di] = water_props(T2d(:),P2d(:),{'H2O,tc-ds633'},'PS94','S14');
v  = null(Nphs);
m0 = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001 0.0001]*100;
for iT = 1:length(T2d(:))
    f = @(m) dmu_fun(m,m_id,g0(:,iT),v,T2d(iT),chrg,rho_w(iT),eps_di(iT));
    m = fsolve(f,m0); % find molalities
    %postprocess
    [gam,a_w] = gam_HKF(T2d(iT),m',chrg(m_id==2),rho_w(iT),eps_di(iT));
    pH(iT) = -log10(m(1)*gam(1)); % pH
    Si(iT) = m(3)+m(4)+2*m(5)+3*m(6)+2*m(7); % Total dissolved Si
    m0 = m;% take previous m as new starting guess
end
figure(1),clf
plot(T-273.15,log10(Si))
% Add Manning 1994 data on top
iP = 1;[data,txt] = xlsread('Manning94_data');T_data = data(:,1);P_data = data(:,2); mSiO2_data = data(:,3);wSiO2_data = data(:,4);
hold on,plot(T_data(P_data == P(iP)/1e8),log10(mSiO2_data(P_data == P(iP)/1e8)),'o','LineWidth',1)
% Function for reactions for fsolve:
function F = dmu_fun(m,m_id,g0,v,T,chrg,rho_w,eps_di)
    [gam,a_w] = gam_HKF(T,m',chrg(m_id==2),rho_w,eps_di); % Activity coefficients for non-ideality from HKF
    a = ones(size(g0)); 
    a(m_id==2) = m'.*gam; % activity of non-pure species
    a(m_id==1) = a_w;     % activity of water
    a(a<0)     = 1e-20;
    for i = 1:length(g0)       
        mu(i) = g0(i) + 8.31446*T*log(a(i)); % definition of chemical potential
    end
    for ieq = 1:size(v,2)
        F(ieq) = (mu*v(:,ieq))/1e3; % chemical reactions with chemical potentials and stoichiometry from nullspace
    end
    F(ieq+1) = m*chrg(m_id==2)'; % charge balance
end