clear,figure(1),clf,colormap jet,addpath ../
T   = 25+273.15;
P   = 1e5;
X = {'Si','Al','K','H','O','e'};
water_id = 'H2O,l,NIST';
phases       = {'MICROCLINE,MAXIMUM,supcrt','MUSCOVITE,supcrt','KAOLINITE,supcrt','GIBBSITE,supcrt','AMORPHOUS-SILICA,supcrt','PYROPHYLLITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
x   = logspace(-7,0,600);
y   = logspace(-7,0,600);
td = init_thermo(phases,X);
p  = props_generate(td);
rho_w     = rho_H2O(T,P,'JN91');
eps_di    = eps_H2O(T,P,rho_w,'JN91');
[g0,v0]   = tl_g0(T,P,td,rho_w,eps_di);
[g0,Nphs] = tl_gibbs_energy(T,P,phases,td,p,g0,v0,rho_w,eps_di);
Nphs = Nphs';
react{1} = {'MICROCLINE,MAXIMUM,supcrt','MUSCOVITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{2} = {'KAOLINITE,supcrt','GIBBSITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{3} = {'MICROCLINE,MAXIMUM,supcrt','KAOLINITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{4} = {'MUSCOVITE,supcrt','KAOLINITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{5} = {'MUSCOVITE,supcrt','GIBBSITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{6} = {'AMORPHOUS-SILICA,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{7} = {'PYROPHYLLITE,supcrt','MICROCLINE,MAXIMUM,supcrt','K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
react{8} = {'PYROPHYLLITE,supcrt','KAOLINITE,supcrt',water_id,'K+,supcrt','H+,supcrt','SiO2,aq,supcrt'};
A_react = zeros(length(react),length(phases));
for i_r = 1:length(react)
    id = [];
    for i = 1:length(react{i_r})
        id(i) = find(strcmp(phases,react{i_r}(i)));
    end
    v = null(Nphs(id,:)');        
    A_react(i_r,id) = v;
end
solid_names = phases(1:6);
iPT = 0;
R   = 8.3144;
% Do reactions:
for ix = 1:length(x)
    for iy = 1:length(y)
        for ir = 1:size(A_react,1)
            a   = [1  1  1  1  1  1 1 1 y(iy) x(ix)];
            Keq = prod(a.^A_react(ir,:));            
            dg0(ir)   = sum(g0'.*A_react(ir,:));
            dG(ix,iy,ir) = dg0(ir) + R*T*log(Keq);            
        end
    end
end
alph = stable_reactions(A_react,dG,solid_names); % this replaces linprog in a way
asm_id = zeros(length(x)*length(y),length(phases));
for i = 1:size(asm_id,1)
    asm_id(i,1:length(find(alph(i,:)>0))) = find(alph(i,:)>0);
end
% Plot phase diagram section
figure(1)
tl_psection(log10(x),log10(1./y),X,asm_id,phases,[0*1e10,0*1e10],8)
xlabel('log(H_4SiO_4)'),ylabel('(log(K+)');
title('Activity diagram');
axis square
set(gca,'FontSize',14)
% print('Fig2d','-depsc','-r600')