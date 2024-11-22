clear,clf,addpath ..\ ..\Solutions\ ..\Utilities\ ..\EOS
run_name = 'ternary_NaCl_KCl_H2O_573K';
load(['linprog_run_' run_name]);       
molm      = molmass_fun(Cname);
solv_tol  = 0.5;
fluid     = phs_name(1);%for feldspar use empty: ''
Cf        = zeros(length(X1_nd(:)),length(Cname));
Cs        = zeros(length(X1_nd(:)),length(Cname));
phs_modes = zeros(length(Tnd(:)),length(phs_name));
for iPT = 1:length(Tnd(:))    
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,g,alph] = postprocess_fun(Tnd(iPT),Pnd(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'CORK','S14');
    asm_id(iPT,1:length(phi)) = pc_id;
    phs_modes(iPT,pc_id) = phi;
    fluid_id =  strcmp(phs_name(pc_id),fluid);
    solid_id = ~strcmp(phs_name(pc_id),fluid);
    if sum(fluid_id)>0
        Cf(iPT,:)     = Cwt(:,fluid_id)*phiw(fluid_id)/sum(Cwt(:,fluid_id)*phiw(fluid_id));  %p_out{pc_id(fluid_id)};
        rhof_tab(iPT) = rho(fluid_id)'*phi(fluid_id)/sum(phi(fluid_id));           
        phif_tab(iPT) = phi(fluid_id);
    end  
    if sum(solid_id)>0
        rhos_tab(iPT) = rho(solid_id)'*phi(solid_id)/sum(phi(solid_id));
        phis_tab(iPT) = sum(phi(solid_id));
        Cs(iPT,:)     = Cwt(:,solid_id)*phiw(solid_id)/sum(Cwt(:,solid_id)*phiw(solid_id));
    end
end
asm_id_all          = zeros(length(X1)*length(X2)*length(T)*length(P),size(asm_id,2));
asm_id_all(id,:)    = asm_id;
asm_id              = reshape(asm_id_all,[length(X1)*length(X2),length(T),size(asm_id,2)]);
Cf_all       = zeros(length(X1)*length(X2)*length(T)*length(P),size(Cf,2));
Cf_all(id,:) = Cf;
Cf           = reshape(Cf_all,[length(X1),length(X2),length(T),size(Cf,2)]);
Cs_all       = zeros(length(X1)*length(X2)*length(T)*length(P),size(Cs,2));
Cs_all(id,:) = Cs;
Cs           = reshape(Cs_all,[length(X1),length(X2),length(T),size(Cs,2)]);
phs_modes_all       = zeros(length(X1)*length(X2)*length(T)*length(P),size(phs_modes,2));
phs_modes_all(id,:) = phs_modes;
solid_id            = find(~strcmp(phs_name,fluid));
solid_names         = phs_name(solid_id);
vol_frac_solids     = phs_modes_all(:,solid_id)./sum(phs_modes_all(:,solid_id)+1e-20,2);
solid_names         = solid_names(sum(vol_frac_solids)>0);
vol_frac_solids     = vol_frac_solids(:,sum(vol_frac_solids)>0);
vol_frac_solids     = reshape(vol_frac_solids,[length(X1)*length(X2),length(T),size(vol_frac_solids,2)]);
phs_modes           = reshape(phs_modes_all,[length(X1)*length(X2),length(T),size(phs_modes,2)]);
[X1_2d,X2_2d]    = ndgrid(X1,X2);
[x2d,y2d]        = cart2bary(X1_2d,X2_2d);
asm_id(X1_2d(:)+X2_2d(:)>1,1,:)=nan;% slows down the tl_psection
figure(1)
tl_psection(x2d,y2d,Cname,squeeze(asm_id),phs_name,1);
figure(2)
phs_modes(X1_2d+X2_2d>1,1,:) = nan;
phs_modes = reshape(phs_modes,length(X1),length(X2),length(phs_name));
vol_frac_solids(X1_2d+X2_2d>1,1,:) = nan;
vol_frac_solids = reshape(vol_frac_solids,length(X1),length(X2),length(solid_names));
subplot(211)
nrow = floor(sqrt(length(solid_names)));
ncol = ceil(length(solid_names)/nrow);
for ip = 1:length(solid_names)
    colormap jet
    subplot(nrow,ncol,ip),pcolor(x2d,y2d,vol_frac_solids(:,:,ip));colorbar,shading flat,
    colorbar,shading flat,axis([0 1/sind(60) 0 1]),axis off,axis square,title(solid_names(ip))
    clim([0,1])
end
figure(3),clf,colormap jet
ncol = floor(sqrt(length(Cname)));
nrow = ceil(length(Cname)/ncol);
for ic = 1:length(Cname)
    subplot(nrow,ncol,ic),pcolor(x2d,y2d,       squeeze(Cf(:,:,1,ic))),title(Cname(ic)),colorbar,shading flat,axis([0 1/sind(60) 0 1]),axis off,axis square
end
figure(4),clf,colormap jet
ncol = floor(sqrt(length(Cname)));
nrow = ceil(length(Cname)/ncol);
for ic = 1:length(Cname)
    subplot(nrow,ncol,ic),pcolor(x2d,y2d,       squeeze(Cs(:,:,1,ic))),title([Cname(ic) ' solid']),colorbar,shading flat,axis([0 1/sind(60) 0 1]),axis off,axis square
end