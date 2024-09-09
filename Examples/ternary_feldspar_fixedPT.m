clear,clf,addpath ..\ ..\Utilities\ ..\Solutions\ ..\EOS
run_name = 'ternary_Feldspar_800C_2';
T = 800 + 273.15;
P = 1e8;
X1 = linspace(0,1,60);
X2 = linspace(0,1,61);
phs_name = {'Feldspar'};
Cname = {'Si','Al','Na','K','Ca','O'};
Nsys0 = [ 3    1    0    1   0    8];
Nsys1 = [ 2    2    0    0   1    8];
Nsys2 = [ 3    1    1    0   0    8];
td = init_thermo(phs_name,Cname,'Metabasite');
td(1).dz(:) = 1/6;
p  = props_generate(td);
[X1_nd,X2_nd,Tnd,Pnd] = ndgrid(X1,X2,T,P);
id = X1_nd(:)+X2_nd(:)<=1;
Tnd = Tnd(id);Pnd = Pnd(id);X1_nd = X1_nd(id);X2_nd = X2_nd(id);
parfor iPT = 1:length(X1_nd(:))
    Nsys = Nsys0*X1_nd(iPT) + Nsys1*X2_nd(iPT) + Nsys2*(1-X1_nd(iPT)-X2_nd(iPT));
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(Tnd(iPT),Pnd(iPT),Nsys,phs_name,p,td);
end
save(['linprog_run_' run_name],'-v7.3');