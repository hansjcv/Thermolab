clear,addpath ../ ../Utilities/ ../Solutions/ ../EOS
T       = linspace(300,600,4) + 273.15;
P       = 0.1e9;
Cname   = {'Al','F','H','O','e'};
Cname_f = {'F'};
Csys    = logspace(-4,-0.2,20);
solids  = {'cor'};
solvent = {'H2O'};
spcs    = {'H+','OH-','HF,aq','F-','Al(OH)2F,aq','AlF+2','Al(F)2+','Al(F)3,aq','Al(F)4-','Al(OH)2(F)2-','HF2-','Al(OH)(F)2,aq','Al(OH)4-','Al(OH)2+','Al(OH)3,aq','AlOH+2','Al+3'};%,
% syms_leg = {''}

%% Select dataset and solvent models
db_species    = 'supcrt';
db_solids     = 'tc-ds55';
db_solvent    = 'tc-ds55';
solvent_model = 'IAPWS';
dielect_model = 'S14';
gam_model     = 3;
%% Preprocess
phase   = [solids,solvent,spcs];
m_id    = [zeros(size(solids)),ones(size(solvent)),ones(size(spcs))*2];
for ip = 1:numel(phase)
    if m_id(ip) == 0
        phase{ip} = [phase{ip},',',db_solids];
    elseif m_id(ip) == 1
        phase{ip} = [phase{ip},',',db_solvent];
    elseif m_id(ip) == 2
        phase{ip} = [phase{ip},',',db_species];
    end    
end
%% Gibbs energy
td             = init_thermo(phase,Cname);
[T2d,P2d]      = ndgrid(T,P);
[g0,Nphs]      = tl_gibbs_energy(T2d(:),P2d(:),phase,td);
[rho_w,eps_di] = water_props(T2d(:),P2d(:),solvent,solvent_model,dielect_model);
chrg           = Nphs(end,:);
%% Change basis
ibasis = [1,2,3,5];
N_phs  = [];
Npc    = Nphs;
for ip    = 1:size(Npc,2)   % ???
    react = null(Npc(:,[ibasis,ip]));%,'r'
    N_phs = [N_phs -react(1:end-1,1)/react(end)];
end
Nphs = N_phs;
Nphs(Nphs<1e-10&Nphs>-1e-10)=0;
%% Reactions
niter          = 200;
max_Res        = 1e-13;
v              = null(Nphs,'rational')';
no_dof         = sum(m_id==2)-size(v,1)-1;
disp(['Degrees of freedom = ' num2str(no_dof)]);
disp_reactions([solids,solvent,spcs],v');
mass = Nphs(4,m_id==2);
act_solids = ones(numel(solids),length(T));
%% Speciation calculations
for iX = 1:length(Csys)
    for iT = 1:length(T)
        m    = 0.1*ones(numel(phase(m_id==2)),1);
        gam  = ones(size(m));
        a_w  = 1;
        for iter = 1:niter
            act    = [act_solids(:,iT);a_w; m.*gam];
            mu     = g0(:,iT)/8.3144/T(iT) + log(act);
            Res    = -[v*mu            ;chrg(m_id==2)*m ; -Csys(iX) + mass*m];
            dRes   =  [v(:,m_id==2)./m';chrg(m_id==2)   ;          mass ];
            del_m  = dRes\Res;% 
            % del_m  = pinv(dRes)*Res;
            del_UR = 1./max(1,-del_m./(1/2.1*m));
            m      = m + del_UR.*del_m;% m      = m + del_m;% m(m<0) = 1e-20;
            z      = chrg(m_id==2)';
            I      = m'*z.^2/2;
            Lgam   = -log10(1+18.01528e-3*sum(m)); % large gamma (for 1 kg of water) conversion factor see Helgeson L gamma page 1293 eq 122, same as eq 2 in Miron 2016, assuming njw mole amount of water-solvent = 1
            a0     = 4.5e-8;
            A      = (1.82483e6)*sqrt(rho_w(iT)*1e-3) ./ (T(iT).*eps_di(iT)).^(3/2);
            B      = (50.2918649e8)*sqrt(rho_w(iT)*1e-3) ./ sqrt(T(iT).*eps_di(iT));
            b_gam = 0.03;
            if gam_model == 1 % Debye-Hückel
                lgam       = -(A*z.^2.*sqrt(I))./(1+a0*ones(size(z)).*B.*sqrt(I));
            elseif gam_model == 2 % Davies
                lgam = -A*z.^2.*(sqrt(I)./(1+sqrt(I))-0.2*I);
            elseif gam_model == 3 % HKF
                lgam    = -(A*z.^2.*sqrt(I))./(1+a0*ones(size(z)).*B.*sqrt(I)) + b_gam*I + Lgam; % probably molar scale see Helgeson L gamma page 1293 eq 122
            end
            gam = 10.^(lgam);
            % Activity of water
            nu        = 2;
            lam       = 1 + a0*B*sqrt(I); % lamda
            sig       = 3./(a0^3*B.^3*sqrt(I.^3)).*(lam-1./lam-2*log(lam));  % sigma coefficient
            phi       = -log(10)*sum(m(z~=0))./sum(m).*(A.*sqrt(I).*sig/3 + Lgam./(18.01528e-3*nu*I) - b_gam*I/2); % Osmotic coefficient
            phi(I==0) = 0; % set osmotic coefficient to 0 for case of no electrolytes
            a_w       = exp(-phi.*sum(m)*18.01528e-3); % Activity of water
            if max(abs(Res(:)))<max_Res,break,end
        end
        iters(iX,iT)   = iter;
        chk(iX,iT)     = max(abs(Res(:)));
        pH(iX,iT)      = -log10(m(strcmp(spcs,'H+')));
        N_tot(iX,iT,:)    = Nphs(:,m_id==2)*m;
        m_all(iX,iT,:) = m;
    end
end
%% Postprocess
max(abs(chk(:)))
[molality,id] = sort(squeeze(m_all(end,1,:)),'descend');
species = spcs(id)';
table(species,molality)
disp(['pH = ' num2str(pH(1,end))])
% figure(1),clf,plot(log10(m_all(:,:,3)),log10(1*N_tot(:,:,1))),axis([-4 0 -5 -1])
figure(1),clf,plot(log10(Csys),log10(2*N_tot(:,:,1))),axis([-4 0 -5 -1])
load zaraisky94_cor_HF_300.mat
hold on
plot(exp_data(:,1),exp_data(:,2),'ok')
load zaraisky94_cor_HF_400.mat
hold on
plot(exp_data(:,1),exp_data(:,2),'vk')
load zaraisky94_cor_HF_500.mat
hold on
plot(exp_data(:,1),exp_data(:,2),'dk')
load zaraisky94_cor_HF_600.mat
hold on
plot(exp_data(:,1),exp_data(:,2),'*k')
legend('300\circC','400\circC','500\circC','600\circC','Zaraisky (1994) at 300 \circC','Zaraisky (1994) at 400 \circC','Zaraisky (1994) at 500 \circC','Zaraisky (1994) at 600 \circC','Location','NorthWest')
xlabel('log molality HF')
ylabel('log total Al')

figure(2),clf,plot(log10(Csys),log10(2*N_tot(:,1,1)),log10(Csys),squeeze(log10(m_all(:,1,:))));axis([-4 0 -5 -1])
load zaraisky94_cor_HF_300.mat
hold on
plot(exp_data(:,1),exp_data(:,2),'ok')
ax = gca;
mylinestyles = ["-"; "--"; "-o"];
ax.LineStyleOrder = mylinestyles;
% load zaraisky94_cor_HF_400.mat
% hold on
% plot(exp_data(:,1),exp_data(:,2),'vk')
% load zaraisky94_cor_HF_500.mat
% hold on
% plot(exp_data(:,1),exp_data(:,2),'dk')
% load zaraisky94_cor_HF_600.mat
% hold on
% plot(exp_data(:,1),exp_data(:,2),'*k')
legend(['total Al',spcs,'Zaraisky (1994) experiments'],'Location','NorthWest')
xlabel('log molality HF')
ylabel('log molality')