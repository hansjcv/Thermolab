clear,figure(1),clf,colormap(jet(256)),addpath ../
% Physics
Lx = 1;                                                        % Model domain length
Dc = 1;                                                        % Diffusion coefficient
% Thermodynamics
T        = 500 + 273.15;                                       % Temperature (K)
P        = 1e8;                                                % Pressure (Pa)
R        = 8.31;
phs_name = {'Orthopyroxene'};                                   % phase name
Cname    = {'Si','Ca','Mg','O'};                           % system components
td       = init_thermo(phs_name,Cname,'solution_models_H18');  % load static thermodynamic data
[g0,V0]  = tl_g0(T,P,td);                                      % get endmember gibbs energies of the solution model
% Numerics
nx   = 200;                                                    % Grid resolution
nt   = 1e6;                                                    % number of time steps
nout = 5000;                                                   % plotting each nout time step
gam  = 2*(Lx/(nx-1))^2*10;                                        % Cahn Hilliard surface energy parameter
delC = 1e-3;                                                   % delta concentration for numerical derivative
% Preprocessing 
dx = Lx/(nx-1);                                                % grid step
dt = min(dx^2/Dc/2.1,dx^4/gam/4.1)/8;                          % time step
x  = 0:dx:Lx;                                                  % x coordinates
% Initialization
c    = 0.5 + 0.01*rand(1,nx);                                  % initial concentration of albite
cini = c;                                                      % store initial concentration for plot
% Processing
for it = 1:nt
    p{1}        = [c; 1-c]';                                   % proportions of endmembers Ab-Or
    p_dC{1}     = [c+delC; 1-(c+delC)]';                       % proportions of endmembers Ab-Or plus delta C
    g_p         = tl_gibbs_energy(T,P,phs_name,td,p_dC,g0,V0); % Gibbs energy + delta G, for numerical derivative
    g           = tl_gibbs_energy(T,P,phs_name,td,p   ,g0,V0); % Gibbs energy
    mu          = (g_p-g)'/delC/R/T;                           % mu1 - mu2; numerical derivative dg/dc, eq. 64
    mu(2:end-1) = mu(2:end-1) - gam*diff(c,2)/dx^2;            % Cahn-Hilliard (1958), eq. 64
    cc          = (c(1:end-1)+c(2:end))/2;                     % Centered concentrations between gridpoints
    qc          = -Dc*(1-cc).*cc.*diff(mu)/dx;                 % Non-Fickean Diffusion flux, eq. 63
    c(2:end-1)  = c(2:end-1) - dt*diff(qc)/dx;                 % Concentration balance with assumptions, eq. 62
    c(1:3)      = mean(c(1:3));                                % no flux left boundary
    c(end-2:end)= mean(c(end-2:end));                          % no flux right boundary
    % Postprocessing
    if mod(it,nout) == 1       
        plot(x,c,x,cini),title([it,mean(c)]),axis([0 1 0 1])   % plot concentration versus distance        
        drawnow
        irec          = (it-1)/nout + 1;
        cevol(:,irec) = c;
        time(irec)    = it*dt;
    end
end
pcolor(log(time),x,cevol),shading flat,colorbar
xlabel('log(time)'),ylabel('x'),title('C(t,x)'),xlim([-4 0.1])