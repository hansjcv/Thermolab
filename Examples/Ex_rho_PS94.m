clear,figure(1),clf,colormap jet,addpath ../ ../EOS/
T = linspace(25,1000,180) + 273.15;
P = linspace(1,2e9,181);
[T2d,P2d] = ndgrid(T,P);
td        = init_thermo({'H2O,tc-ds633','CO2,tc-ds633'});
[g0,v0]   = tl_g0(T2d(:),P2d(:),td);
molmW     = [18.0153e-3,44.0095e-3];
rho_kgm3(:,1)  = molmW(1)./v0{1}*1e5;
rho_kgm3(:,2)  = molmW(2)./v0{2}*1e5;
subplot(121),contourf(T2d-273.15,P2d/1e9,reshape(rho_kgm3(:,1),length(T),length(P)),30),colorbar,shading flat
title('\rho_{H_2O} (kg/m^3)'),xlabel('T(\circC)'),ylabel('P(GPa)'),caxis([0 1650])
subplot(122),contourf(T2d-273.15,P2d/1e9,reshape(rho_kgm3(:,2),length(T),length(P)),30),colorbar,shading flat
title('\rho_{CO_2} (kg/m^3)'),xlabel('T(\circC)'),ylabel('P(GPa)'),caxis([0 1650])
