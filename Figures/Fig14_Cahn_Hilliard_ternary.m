clear,clf,addpath ../
load run_ternary_25_10_2021
imagesc(c1),title('ab','FontSize',18),colorbar,caxis([0 1]),set(gca,'FontSize',18),axis square,axis off
% print('Fig14c_CH_ternary_2D','-r600','-dpng')
imagesc(c2),title('an','FontSize',18),colorbar,caxis([0 1]),set(gca,'FontSize',18),axis square,axis off
% print('Fig14d_CH_ternary_2D','-r600','-dpng')
imagesc(1-c1-c2),title('or','FontSize',18),colorbar,caxis([0 1]),set(gca,'FontSize',18),axis square,axis off
% print('Fig14e_CH_ternary_2D','-r600','-dpng')
[c1_2d,c2_2d] = cart2bary(c1_plot,c2_plot);
[c1_p,c2_p] = cart2bary(c1,c2);
p_labels = {'ab','an','or'};
p1_end   = [1    0    0];
p2_end   = [0    1    0];
[c1_end,c2_end] = cart2bary(p1_end,p2_end);
[c1_sys,c2_sys] = cart2bary(psys(1),psys(2));


run_name = 'ternary_fsp_2021_10_13';
load(['linprog_run_' run_name]);
nstab = zeros(length(alph_all),1);
[x2d,y2d] = cart2bary(c12d,c22d);
solv_tol = 0.075;
for iPT = 1:length(alph_all)    
    if ~isnan(alph_all{iPT})
        [N_mol,Npc,p,pc_id] = cluster_p(alph_all{iPT},Npc_all{iPT},p_all{iPT},pc_id_all{iPT},solv_tol,phs_name);
        nstab(iPT) = length(pc_id);
        pstab = p_all{iPT}{1}(alph_all{iPT}>0,:);
    else
        nstab(iPT) = nan;
    end       
end
[p1,p2] = cart2bary(c12d,c22d);
contourf(c1_2d,c2_2d ...
    ,reshape(g_plot-g_mech,size(c1_plot))),axis equal,axis off
hold on,
contour(p1,p2,reshape(nstab,length(c1),length(c2)),[1.5,1.5],'y--','LineWidth',3)
axis([0 1.2 0 1])
plot(c1_p,c2_p,'yo'),plot(p1eq,p2eq,'rd','MarkerSize',10,'LineWidth',4)

plot(c1_sys,c2_sys,'yo','MarkerSize',10,'LineWidth',4)
text(c1_end+[0.05 -0.03 -0.1],c2_end+[-0.05 0.05 -0.05],p_labels,'FontSize',18)%,hold off
colorbar
colormap jet
set(gcf,'Renderer','ZBuffer')
% print('Fig14f_CH_ternary_2D','-r600','-dpng')