clear,figure(1),clf,colormap jet,addpath ../
% Physics
Lx = 1;
Ly = 1;
Dc = 1;
% Thermodynamics
T        = 500 + 273.15;
P        = 1e8;
phs_name = {'Feldspar(C1)'};
Cname    = {'Si','Al','Ca','Na','K','O'};
% Numerics
nx   = 150;
ny   = 150;
nt   = 1000000;
nout = 100;
delC = 1e-3;
% Preprocessing
dx      = Lx/(nx-1);
dy      = Ly/(ny-1);
gam     = 4e-3*dx^2;
dt      = min(min(dx^2,dy^2)/Dc/2.1)/1;
x       = 0:dx:Lx;
y       = 0:dy:Ly;
td      = init_thermo(phs_name,Cname,'solution_models_HP98'); % load static thermodynamic data
[g0,V0] = tl_g0(T,P,td); % get endmember gibbs energies of the solution model
% Initialization
c1      = 0.15 + 0.25*rand(nx,ny);
c2      = 0.15 + 0.25*rand(nx,ny);
qcx1     = zeros(nx+1,ny);
qcx2     = zeros(nx+1,ny);
qcy1     = zeros(nx,ny+1);
qcy2     = zeros(nx,ny+1);
% Use only to plot Gibbs mixing curve
[c12_plot,c3_plot]  = ndgrid(linspace(0,1,100));
c1_plot   =    c12_plot .*(1-c3_plot);
c2_plot   = (1-c12_plot).*(1-c3_plot);
p_plot{1} = [c1_plot(:),c2_plot(:),c3_plot(:)];
g_plot    = tl_gibbs_energy(T,P,phs_name,td,p_plot,g0,V0);
gscale    = 1/abs(max(g_plot));
g_mech    = p_plot{1}*g0{1}';
contourf(c1_plot+c2_plot/2,c2_plot*sqrt(3)/2,reshape(g_plot-g_mech,size(c1_plot))),axis equal,axis off
surf(c1_plot+c2_plot/2,c2_plot*sqrt(3)/2,reshape(g_plot-g_mech,size(c1_plot))),shading flat,light,lighting phong
% Processing
for it = 1:nt
    c1(c1>1-2*delC) = 1-2*delC; c1(c1<delC) = delC;
    c2(c2>1-2*delC) = 1-2*delC; c2(c2<delC) = delC;    
    p{1}            = [c1(:),c2(:), 1-c1(:)-c2(:)]; % proportions of endmembers Ab-An-San
    p_dC1{1}        = [c1(:)+delC, c2(:)     , 1-(c1(:)+c2(:)+delC)];
    p_dC2{1}        = [c1(:)     , c2(:)+delC, 1-(c1(:)+c2(:)+delC)];
    g_p1            = tl_gibbs_energy(T,P,phs_name,td,p_dC1,g0,V0);
    g_p2            = tl_gibbs_energy(T,P,phs_name,td,p_dC2,g0,V0);
    g               = tl_gibbs_energy(T,P,phs_name,td,p    ,g0,V0);
    mu1             = reshape(gscale*(g_p1-g)/delC,nx,ny); % mu1 - mu3; numerical derivative dg/dc1
    mu2             = reshape(gscale*(g_p2-g)/delC,nx,ny); % mu2 - mu3; numerical derivative dg/dc2
    mu1(2:end-1,:)   = mu1(2:end-1,:) - gam*diff(c1   + c2/1,2,1)/dx^2; % Cahn-Hilliard
    mu2(2:end-1,:)   = mu2(2:end-1,:) - gam*diff(c1/1 + c2  ,2,1)/dx^2; % Cahn-Hilliard
    mu1(:,2:end-1)   = mu1(:,2:end-1) - gam*diff(c1   + c2/1,2,2)/dx^2; % Cahn-Hilliard
    mu2(:,2:end-1)   = mu2(:,2:end-1) - gam*diff(c1/1 + c2  ,2,2)/dx^2; % Cahn-Hilliard
    mu1([1 (end)],:) = mu1([2 (end-1)],:);
    mu2([1 (end)],:) = mu2([2 (end-1)],:);
    mu1(:,[1 (end)]) = mu1(:,[2 (end-1)]);
    mu2(:,[1 (end)]) = mu2(:,[2 (end-1)]);
    qcx1(2:end-1,:)   = -Dc*( diff(mu1,1,1)/dx   - diff(mu2,1,1)/dx/2); % Non-Fickean Diffusion flux
    qcx2(2:end-1,:)   = -Dc*(-diff(mu1,1,1)/dx/2 + diff(mu2,1,1)/dx  ); % Non-Fickean Diffusion flux
    qcy1(:,2:end-1)   = -Dc*( diff(mu1,1,2)/dy   - diff(mu2,1,2)/dy/2); % Non-Fickean Diffusion flux
    qcy2(:,2:end-1)   = -Dc*(-diff(mu1,1,2)/dy/2 + diff(mu2,1,2)/dy  ); % Non-Fickean Diffusion flux
    c1             = c1 - dt*diff(qcx1,1,1)/dx; 
    c2             = c2 - dt*diff(qcx2,1,1)/dx;
    c1             = c1 - dt*diff(qcy1,1,2)/dy; 
    c2             = c2 - dt*diff(qcy2,1,2)/dy;
    % Postprocessing
    if mod(it,nout) == 1
        g_m    = p{1}*g0{1}';
        subplot(221),imagesc(c1),axis square
        subplot(222),imagesc(c2),axis square
        subplot(223),imagesc(1-c1-c2),axis square
        subplot(224),contourf(c1_plot+c2_plot/2,c2_plot*sqrt(3)/2 ...
            ,reshape(g_plot-g_mech,size(c1_plot))),axis equal,axis off
        hold on,plot(c1+c2/2,c2*sqrt(3)/2,'wo'),hold off,title(it)
        drawnow
    end
end