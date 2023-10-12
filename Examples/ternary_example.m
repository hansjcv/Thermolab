clear,clf,addpath ..\ ..\Utilities\ ..\Solutions\
run_name = 'ternary_NaCl_KCl_H2O_573K';
T = 300 + 273.15;
P = 0.5e8;
X1 = linspace(0,1,60);
X2 = linspace(0,1,61);
phs_name = {'Salt(L)','Salt'};
Cname = {'Na','K','Cl','H','O'};
Nsys0 = [ 1    0   1    0    0];
Nsys1 = [ 0    0   0    2    1];
Nsys2 = [ 0    1   1    0    0];
td = init_thermo(phs_name,Cname,'solution_models_HP98');
td(1).dz(:) = 1/100;
p  = props_generate(td);
[X1_nd,X2_nd,Tnd,Pnd] = ndgrid(X1,X2,T,P);
id = X1_nd(:)+X2_nd(:)<=1;
Tnd = Tnd(id);Pnd = Pnd(id);X1_nd = X1_nd(id);X2_nd = X2_nd(id);
parfor iPT = 1:length(X1_nd(:))
    Nsys = Nsys0*X1_nd(iPT) + Nsys1*X2_nd(iPT) + Nsys2*(1-X1_nd(iPT)-X2_nd(iPT));
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(Tnd(iPT),Pnd(iPT),Nsys,phs_name,p,td);
end
save(['linprog_run_' run_name],'-v7.3');