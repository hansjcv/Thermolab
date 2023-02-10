clear,clf
runname = 'K32_2023_01_26_H18';
load(['linprog_run_' runname]);                                          % load linprog run data
molm = molmass_fun(Cname);
solv_tol = 2;
fluid = 'Fluid-H2O';
phs_modes = zeros(length(T2d(:)),length(phs_name));
for iPT = 1:length(T2d(:))    
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,mu2] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'CORK','S14');
    phs_modes(iPT,pc_id) = phi;
    asm_id(iPT,1:length(phi)) = pc_id;
end
figure(1),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name,0,[0 0],6);
figure(2),
nrow = floor(sqrt(length(phs_name)));
ncol = ceil(length(phs_name)/nrow);
phs_modes = reshape(phs_modes,length(T),length(P),length(phs_name));
for ip = 1:length(phs_name)
    colormap jet
    subplot(nrow,ncol,ip),pcolor(T2d-273.15,P2d/1e9,phs_modes(:,:,ip));colorbar,shading flat,title(phs_name(ip))
end