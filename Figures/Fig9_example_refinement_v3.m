clear,close all,addpath ../
xsys = 0.38;
g0A = [0    0.1];
g0B = [0.15 0  ];
xA_in = [0
    1];
xB_in = [0
    1];
nref = 2;
for i_ref = 1:nref
    % Grid
    xA = linspace(xA_in(1),xA_in(2),10);dxA = xA(2)-xA(1);
    xB = linspace(xB_in(1),xB_in(2),10);dxB = xB(2)-xB(1);
    % Gibbs energy
    gA = g0A(1)*xA + g0A(2)*(1-xA) + xA.*log(xA+double(xA==0)) + (1-xA).*log(1-xA+double(xA==1));
    gB = g0B(1)*xB + g0B(2)*(1-xB) + xB.*log(xB+double(xB==0)) + (1-xB).*log(1-xB+double(xB==1));
    % Minimization
    Aeq  = [xA xB; 1-xA 1-xB];
    beq  = [xsys;1-xsys];
    g    = [gA gB];
    LB   = zeros(1,size(g,2));
    alph = linprog(g,[],[],Aeq,beq,LB);
    % Postprocess
    xstable(:,i_ref) = Aeq(1,alph>0);
    gstable = g(alph>0);
    % Refinement
    xA_in = [xstable(1,i_ref)-dxA;xstable(1,i_ref)+dxA];
    xB_in = [xstable(2,i_ref)-dxA;xstable(2,i_ref)+dxB];
    xA_in(xA_in<0) = 0;xA_in(xA_in>1) = 1;
    xB_in(xB_in<0) = 0;xB_in(xB_in>1) = 1;
    % Plotting
%     subplot(2,1,i_ref),
    figure(i_ref)
    plot(xA,gA,'bo-',xB,gB,'r*-','LineWidth',1),hold on,
    plot(xstable(:,i_ref),gstable,'go-','LineWidth',2,'MarkerSize',10)   
    xlabel('X'),ylabel('G'),title(['number of refinements = ' num2str(i_ref)])
    legend({'phase A','phase B','tangent'} )   
    axis square
%     print(['Fig8_' num2str(i_ref)],'-depsc','-r600')
end
xstable_ref = xstable(:,end)'
load brute_xstable
abs(xstable-xstable_ref)