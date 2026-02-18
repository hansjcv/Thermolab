clear,clf,addpath ../ ../EOS
T        = linspace(25,1400,30)+273.15;
P        = linspace(0.1,4,31)*1e9;
Tr = 25 + 273.15;
phases       = {'atg,tc-ds62','br,tc-ds62','fo,tc-ds62','H2O,tc-ds62','diam,tc-ds62','gph,tc-ds62','coe,tc-ds62','q,tc-ds62','an,tc-ds62','gr,tc-ds62','ky,tc-ds62'};
i_react(1,:) = [1      1    1    1     0      0     0    0    0    0    0 0 0 0 0];symcol{1}='b';leg_name{1} = 'atg + br = fo + H2O';
i_react(2,:) = [0      0    0    0     1      1     0    0    0    0    0 0 0 0 0];symcol{2}='g';leg_name{2} = 'diam = grph';
i_react(3,:) = [0      0    0    0     0      0     1    1    0    0    0 0 0 0 0];symcol{3}='r';leg_name{3} = 'coe = q';
i_react(4,:) = [0      0    0    0     0      0     0    1    1    1    1 0 0 0 0];symcol{4}='c';leg_name{4} = 'an = ky + gr + q';
[T2d,P2d] = ndgrid(T,P);
td = init_thermo(phases);
[G,Nphs] = tl_gibbs_energy(T2d(:),P2d(:),phases,td);
for i_r = 1:size(i_react,1)
    v   = null(Nphs(:,i_react(i_r,:)==1),'r');
    iph = find(i_react(i_r,:)==1);
    dG(:,i_r)  = G(iph,:)'*v;
    contour(T2d-273.15,P2d/1e9,reshape(dG(:,i_r),length(T),length(P)),[0,0],symcol{i_r},'LineWidth',1)%,axis([200 800 0.5 4]),xlabel('T\circC'),ylabel('P(GPa)')
    hold on
end
legend(leg_name,'Location','NorthWest')
axis square
set(gca,'FontSize',10)