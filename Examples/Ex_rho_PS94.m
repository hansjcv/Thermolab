clear,clf,addpath ../ ../EOS/
T = linspace(25,600,180) + 273.15;
P = linspace(1,0.01e8,181);
[T2d,P2d] = ndgrid(T,P);
td        = init_thermo({'H2O,tc-ds633'});
[g0,v0]   = tl_g0(T2d(:),P2d(:),td);
molmW     = 18.0153e-3;
rho_kgm3  = molmW./v0{1}*1e5;
figure(1),colormap jet
contourf(T2d-273.15,P2d/1e9,reshape(rho_kgm3,length(T),length(P)),30),colorbar,shading flat