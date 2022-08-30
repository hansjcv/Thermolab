clear,close all,addpath ../
iXp = 100;
font_size = 14;
load Fig8a_ol_diagram_dz15.mat
% Figure a
figure(1),clf,colormap jet
pcolor(Xsys,T-273.15,melt_frac'),colorbar,shading flat,title([{'melt amount'};{'n = 15'}]),ylabel('T(\circC)'),axis square
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
set(gca,'Fontsize',font_size)
% print('Fig_8a','-depsc','-r600')
% Figure b
figure(2),clf,colormap jet
pcolor(Xsys,T-273.15,fo_melt'),colorbar,shading flat,title([{'melt composition'},{'n = 15'}]),xlabel('X_{sys}'),axis square
set(gca,'Fontsize',font_size)
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
% print('Fig_8b','-depsc','-r600')
% Figure h
figure(3),clf,colormap jet
plot(T-273.15,fo_melt(iXp,:),'.-b','LineWidth',1)%,title('n = 15')
xlabel('T(\circC)'),ylabel('X_{fo}')
axis square
set(gca,'Fontsize',font_size)
% Figure g
figure(6),clf,colormap jet
plot(T-273.15,melt_frac(iXp,:),'.-b','LineWidth',1)%,title('n = 15')
xlabel('T(\circC)'),ylabel('mol abundance')
axis square
set(gca,'Fontsize',font_size)
load Fig8c_ol_diagram_dz50.mat
% Figure c
figure(4),clf,colormap jet
pcolor(Xsys,T-273.15,melt_frac'),colorbar,shading flat,title('n = 50')%,xlabel('X_{sys}'),
ylabel('T(\circC)'),
axis square
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
set(gca,'Fontsize',font_size)
% print('Fig_8c','-depsc','-r600')
% Figure d
figure(5),clf,colormap jet
pcolor(Xsys,T-273.15,fo_melt'),colorbar,shading flat,title('n = 50')%,xlabel('X_{sys}'),ylabel('T(\circC)'),
axis square
set(gca,'Fontsize',font_size)
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
% print('Fig_8d','-depsc','-r600')
figure(3),hold on
plot(T-273.15,fo_melt(iXp,:),'.-r','LineWidth',1)%,title('melt comp')
xlabel('T(\circC)'),ylabel('X_{fo}')
axis square
set(gca,'Fontsize',font_size)
figure(6),hold on
plot(T-273.15,melt_frac(iXp,:),'.-r','LineWidth',1)%,title('melt fraction')
xlabel('T(\circC)'),ylabel('mol abundance')
axis square
set(gca,'Fontsize',font_size)
iXp = 100;
% load Fig8e_ol_diagram_dz100.mat
% Figure e
figure(7),clf,colormap jet
pcolor(Xsys,T-273.15,melt_frac'),colorbar,shading flat,title('n = 100'),
xlabel('X_{sys}'),
ylabel('T(\circC)'),axis square
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
set(gca,'Fontsize',font_size)
% print('Fig_8e','-depsc','-r600')
% Figure f
figure(8),clf,colormap jet
pcolor(Xsys,T-273.15,fo_melt'),colorbar,shading flat,title('n = 100')%,ylabel('T(\circC)')
xlabel('X_{sys}'),axis square
set(gca,'Fontsize',font_size)
hold on
plot(ones(size(T))*Xsys(iXp),T-273.15,'--w','LineWidth',1)
% print('Fig_8f','-depsc','-r600')
figure(3),hold on
plot(T-273.15,fo_melt(iXp,:),'.-g','LineWidth',1)%,title('melt comp')
xlabel('T(\circC)'),ylabel('X_{fo}')
axis square
set(gca,'Fontsize',font_size)
legend('n = 15','n = 50','n = 100')
% print('Fig_8g','-depsc','-r600')
figure(6),hold on
plot(T-273.15,melt_frac(iXp,:),'.-g','LineWidth',1)%,title('melt fraction')
xlabel('T(\circC)'),ylabel('mol abundance')
axis square
set(gca,'Fontsize',font_size)
legend('n = 15','n = 50','n = 100')
% print('Fig_8h','-depsc','-r600')