clear,addpath ..\ ..\Utilities\ ..\Solutions\
run_name = 'ternary_example_Tpz';
T = 400 + 273.15;
P = 1e6;
X1 = linspace(0,1,30);
X2 = linspace(0,1,31);
phs_name = {'Fluid-KCl-NaCl-HF-H2O','Topaz','mic,tc-ds55','mu,tc-ds55','ab,tc-ds55','and,tc-ds55','q,tc-ds55','syv,tc-ds55','hlt,tc-ds55','cor,tc-ds55'};
Cname = {'Si','Al','Na','K','F','Cl','H','O'};
td = init_thermo(phs_name,Cname,'solution_models_HP98_Tpz');
p  = props_generate(td);
Nsys0 = [3 1 1 0 0 0 0 8]; % Ab;
Nsys1 = [3 1 0 1 0 0 0 8]; % Kfsp
Nsys2 = [1 2 0 0 2 0 0 4]; % Tpz
Nq    = [1 0 0 0 0 0 0 2]; % q;
Nw    = [0 0 0 0 0 0 2 1]*0.5; % H2O;
NCl   = [0 0 1 1 0 2 0 0]*0.01;
[X1_nd,X2_nd,Tnd,Pnd] = ndgrid(X1,X2,T,P);
id = X1_nd(:)+X2_nd(:)<=1;
Tnd = Tnd(id);Pnd = Pnd(id);X1_nd = X1_nd(id);X2_nd = X2_nd(id);
parfor iPT = 1:length(X1_nd(:))
    Nsys = Nsys0*X1_nd(iPT) + Nsys1*X2_nd(iPT) + Nsys2*(1-X1_nd(iPT)-X2_nd(iPT)) + Nq + Nw + NCl;
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(Tnd(iPT),Pnd(iPT),Nsys,phs_name,p,td);
end
save(['linprog_run_' run_name],'-v7.3');