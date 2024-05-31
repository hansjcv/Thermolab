function plot_phase_comp(run_name,phs_id,ps_type,solv_tol)
load(['linprog_run_' run_name]);
molm = molmass_fun(Cname);
phs_modes = zeros(length(T2d(:)),length(phs_name));
np = length(td(phs_id).p_name);
phs_comp  = zeros(length(T2d(:)),np);
for iPT = 1:length(T2d(:))
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,mu2] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'CORK','S14');
    asm_id(iPT,1:length(phi)) = pc_id;
    phs_modes(iPT,pc_id)      = phi;
    if sum(pc_id==phs_id)
        phs_comp(iPT,:)       = p_out{phs_id};
    end
end
switch ps_type
    case 'P-T'
        phs_comp = reshape(phs_comp,length(T),length(P),np);
        nrow          = fix(sqrt(np));% Make number of rows for subplot
        ncol          = ceil(np/nrow);% Make number of columns for subplot
        for ipl = 1:np
            subplot(nrow,ncol,ipl),contourf(T2d-273.15,P2d/1e9,phs_comp(:,:,ipl)*100),colorbar
            title(td(phs_id).p_name(ipl)),xlabel('T(\circC)'),ylabel('P(GPa)')
        end
    case 'T-X'
        phs_comp = reshape(phs_comp,length(X),length(T),np);
        nrow          = fix(sqrt(np));% Make number of rows for subplot
        ncol          = ceil(np/nrow);% Make number of columns for subplot
        for ipl = 1:np
            subplot(nrow,ncol,ipl),contourf(X2d,T2d-273.15,phs_comp(:,:,ipl)*100),colorbar
            title(td(phs_id).p_name(ipl)),ylabel('T(\circC)'),xlabel('X')
        end
end