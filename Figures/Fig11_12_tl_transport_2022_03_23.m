clear,figure(1),clf,colormap(jet),addpath ../
% Physics
load lookup_soapstone_2022_02_14_n1500_nz15_full
mu_tab   =   mu_tab/max(abs(  mu_tab));
rhos_tab = rhos_tab/max(abs(rhos_tab));
rhof_tab = rhof_tab/max(abs(rhof_tab));
%independent
Lx     = 1;                                                               % Length of model domain (m)
Dc     = 1;                                                               % Diffusion coefficient (m^2/s)
P1     = 1;                                                               % Boundary fluid pressure (Pa)
%nondim
phi0   = 0.2;                                                             % Background porosity (fluid volume fraction)
C0     = 0.004;                                                           % Background CO2 concentration in solid (wt%)
C1     = 0.04715;                                                         % Csolid soapstone
npow   = 3;                                                               % Non-linearity in permeability-porosity relation
%dependent
k_muf0 = 1*Dc/P1;                                                         % Background permeability/fluid viscosity (m^2/(Pa*s))
t_tot  = 1e+1*Lx^2/Dc;                                                    % Total time (s)
% Numerics
nx     = 200;                                                             % Number of grid points in x direction
nt     = 1e9;                                                             % Number of time steps
niter  = 1e5;
nout   = 10000;                                                           % Maximum number of Pf iterations
tol    = 1e-6;                                                            % Exit tolerance from Pf iterations
% preprocessing
dx     = Lx/(nx-1);                                                       % Grid step in x
x      = 0:dx:Lx;                                                         % X-coordinates of the nodes
% Initialisation;
phi    = x*0 + phi0;                                                      % Initial porosity
Pf     = P1*(1-x/Lx);                                                     % Initial fluid pressure
Cs     = x*0 + C0;                                                        % Initial solid concentration
Cs(1)  = C1;                                                              % Left solid CO2 concentration boundary
time   = 0;
% Processing
for it = 1:nt   % Time loop
    mu      = interp1(Cs_tab,mu_tab,Cs);                                  % Local equilibrium mu CO2 in fluid from lookup table
    Cf      = interp1(Cs_tab,Cf_tab,Cs);                                  % Local equilibrium CO2 in fluid from lookup table
    Cs_im   = interp1(Cs_tab,Mg_tab,Cs);                                  % Local equilibrium Mg  in solid from lookup table
    rhos    = interp1(Cs_tab,rhos_tab,Cs);                                % Local equilibrium density in solid from lookup table
    rhof    = interp1(Cs_tab,rhof_tab,Cs);                                % Local equilibrium density in fluid from lookup table
    if it == 1,phi_srho_sCs_im0 = (1 - phi).*rhos.*Cs_im;end              % Store initial values
    phi     = 1 - phi_srho_sCs_im0./(Cs_im.*rhos);                        % mass balance of immobile species in solid (eq.7)
    Ctot    = Cf.*rhof.*phi + Cs.*rhos.*(1-phi);                          % Shorthand notation (eq. 2)
    rhotot  =     rhof.*phi +     rhos.*(1-phi);                          % Shorthand notation (eq. 6)
    if it == 1,rhotot0         = rhotot;end                               % Store initial values
    drhotot = rhotot(2:end-1) - rhotot0(2:end-1);                         % Shorthand notation
    % averaging
    rhofc   = 0.5*(rhof(1:end-1)              + rhof(2:end)           );  % rhof values between the nodes
    rhofCfc = 0.5*(rhof(1:end-1).*Cf(1:end-1) + rhof(2:end).*Cf(2:end));  % rhof*Cf values between the nodes
    phic    = 0.5*( phi(1:end-1)              +  phi(2:end)           );  % phi values between the nodes
    % transport properties
    perm    = k_muf0*phic.^npow;                                          % Permeability
    Deff    = Dc*rhofCfc.*phic;                                           % Effective diffusivity    
    % Updates
    for iter = 1:niter        
        dt_dif       = 0.45*dx^2/Dc;                                      % maximum time step of diffusion
        dt_adv       = 0.45*dx/max(abs(rhofc.*perm.*diff(Pf)/dx));        % maximum time step of advection
        dt           = min([dt_adv,dt_dif, t_tot-time]);                  % maximum time step for numerical stability        
        qD           = -perm.*diff(Pf)/dx;                                % Darcy flux (eq. 3)
        Pres         = -drhotot/dt - diff(rhofc.*qD)/dx;                  % Residual of the total mass balance (eq. 5)
        Pf(2:end-1)  = Pf(2:end-1) + dx^2/max(perm.*rhofc)/4*Pres;        % Pf iterative update
        if max(abs(Pres)) < tol,break,end                                 % Exit criteria from Pf iterations
    end
    iters(it)        = iter;                                              % store for plotting
    qC               = -Deff.*diff(mu)/dx;                                % Diffusion flux in fluid (first term in rhs eq. 3)
    Ctot(2:end-1)    = Ctot(2:end-1) - dt*(diff(qC + rhofCfc.*qD)/dx);    % mass concentration balance of mobile species (eq.1)
    Cs               = (Ctot - rhof.*Cf.*phi)./rhos./(1-phi);             % calculate new solid concentration (from eq. 2)
    time             = time + dt;
    rhotot0          = rhotot;                                            % for the next time step
    if mod(it,nout)==1 || time == t_tot  %Postprocessing                  % plot every nout time step
        phs_modes   = interp1(Cs_tab,vol_frac_solids,Cs);                 % find local equilibrium stable phase volume fraction        
        area(x,phs_modes,'FaceColor','flat'),axis tight
        legend(solid_names),title([mean(iters)/nx,iter/nx])        
        xlabel('x')
        drawnow
    end
    if time >= t_tot,break,end
end