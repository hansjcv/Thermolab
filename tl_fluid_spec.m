function [g_fl,N_fl,m_all,chk,logK,v] = tl_fluid_spec(T,B,solvent,spcs,indep_id,g0,Nphs,rho_w,eps_di)
% clear,addpath ../ ../Utilities/ ../Solutions/ ../Figures/ ../EOS
% T     = 320 + 273.15;
% P     = 1e8;
% solvent = {'H2O,tc-ds633'};
% spcs    = {'H+,Miron','OH-,Miron','HF,aq,supcrt','Cl-,Miron','HCl0,Miron','K+,Miron','F-,supcrt','KCl0,Miron','KOH0,Miron'};
phase   = [solvent,spcs];
m_id    = [ones(size(solvent)),ones(size(spcs))*2];
% Cname   = {'K','F','Cl','H','O','e'};
% Cl_B    = linspace(0.0001,0.2,100);
% F_B     = 0.1;
% K_B     = 0.1;
Mw    = 18.01528e-3;
niter = 200;
eps_iter = 1e-12;
m0_guess = 1e-12;
% td    = init_thermo(phase,Cname);
% [T2d,P2d]      = ndgrid(T,P);
% [g0,Nphs]      = tl_gibbs_energy(T2d(:),P2d(:),phase,td);
% [rho_w,eps_di] = water_props(T2d(:),P2d(:),{'H2O,tc-ds633'},'PS94','S14');
% indep_id = 1:3;
v    = null(Nphs,'rational')';
chrg = Nphs(end,:);
m0   = m0_guess*ones(size(phase(m_id==2)));
logK    = log10(exp(-v*g0/8.3144/T));
mass    = Nphs(indep_id,m_id==2);
for iX = 1:size(B,2)
    % B     = [K_B,F_B,Cl_B(iX)]';
    m    = m0_guess*ones(size(phase(m_id==2)));
    gam  = ones(size(m));
    a_w  = 1;
    for  iter = 1:niter
        dK   = v(:,m_id==2)./m/log(10);
        s    = chrg(m_id==2);
        C    = [s;mass;dK];
        zm   = s*m';
        if ~isempty(mass)
            bm   = -B(:,iX) + mass*m';
        else
            bm = [];
        end
        a    = [a_w, m.*gam];
        for i = 1:length(a)
            log_a(i) = log10(a(i)); % activity
        end
        for ieq = 1:size(v,1)
            km(ieq) = log_a*v(ieq,m_id~=0)' - logK(ieq);% equilibrium constant
        end
        [gam,a_w] = gam_HKF(T,m',chrg(m_id==2),rho_w,eps_di);
        gam       = gam';
        Y         = -[zm;bm;km'];
        deltaa    = pinv(C)*Y;
        m         = m + deltaa';
        m(m<0)    = 1e-20;
        test = abs((m-m0)./m0).*100;
        m0 = m;
        if max(test)<eps_iter,break,end
    end
    m_all(iX,:) = m;
    chk(iX)     = max(abs(km));
    s           = m./(sum(m)+1/Mw);
    sw          = 1-sum(s);
    N_fl(:,iX)  = Nphs*[sw;s'];
    g_fl(iX) = sum((g0(m_id==2)' + 8.3144*T*log(a(2:end))).*s) + (g0(m_id==1) + 8.3144*T*log(a_w))*sw; % somehow works as well, is it equivalent to the first one? Yes, subsitute z in the first equation and the log(Nw) cancels
    % pH(iX) = -log10(m_all(iX,1));
end
% disp(num2str(max(abs(chk))))
% figure(1),
% subplot(311),semilogy(pH,m_all(:,2:end)),legend(spcs(2:end))
% subplot(312),semilogy(pH,m_all(:,[3,4,6,7])),legend(spcs([3,4,6,7]))
% subplot(313),semilogy(pH,m_all(:,[5,9])),legend(spcs([5,9]))
% figure(2),plot(solid_modes(:,2)*Nphs_s(4,2),F)