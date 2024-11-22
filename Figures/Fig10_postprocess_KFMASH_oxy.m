clear,addpath ../ ../Utilities/ ../EOS
runname = 'KFMASH_2022_03_06_10x11_oxy';
load(['linprog_run_' runname]);                                          % load linprog run data
molm = molmass_oxy(Cname);
solv_tol = 1;
fluid = 'H2O,tc-ds55';
for iPT = 1:length(T2d(:))    
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,mu2] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'CORK','S14');    
    asm_id(iPT,1:length(phi)) = pc_id;
end
figure(2),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name,0,[0 0],10);