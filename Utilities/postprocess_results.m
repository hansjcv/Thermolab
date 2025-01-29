function [rhos,rhof,cwt_solid,cwt_fluid,phs_modes,Nsfu,asm_id] = postprocess_results(runname,fluid,ndim,ps_type,solv_tol2)
load(['linprog_run_' runname]);                                          % load linprog run data
if nargin>4
    solv_tol = solv_tol2;
end
% Cname    = {'SiO2','Al2O3' , 'CaO' ,'Cr2O3','H2O'};
% Cname    = {'Si','Al' , 'Ca' ,'H' 'O'};
n_ox = zeros(size(Cname));
for i = 1:length(Cname)
    ox_string = Cname{i};
    digit_id  = isstrprop(ox_string,'digit');
    ox_letter = ox_string(~digit_id);
    upper_id  = find(isstrprop(ox_letter,'upper'));
    if numel(upper_id)>1
        cat_name{i} = ox_letter(upper_id(1):upper_id(2)-1);
    else
        cat_name{i} = ox_letter(upper_id(1):end);
        n_cat(i) = 1;
        continue
    end
    nele_id     = find(digit_id);   
    if sum(digit_id)==0 % in case no numbers are found
        n_cat(i)    = 1;
        n_ox(i)     = 1;
    elseif numel(nele_id)==2
        n_cat(i)    = eval(ox_string(nele_id(1)));
        n_ox(i)     = eval(ox_string(nele_id(2)));
    elseif numel(nele_id)==1 && strfind(Cname{i},'O')~=numel(ox_string) % now there are more than one O's
        n_cat(i)    = 1;
        n_ox(i)     = eval(ox_string(nele_id(1)));
    elseif numel(nele_id)==1 && strfind(Cname{i},'O')==numel(ox_string) % now there are more than one O's
        n_cat(i)    = eval(ox_string(nele_id(1)));
        n_ox(i)     = 1;
    end
end
molm = (n_cat.*molmass_fun(cat_name)' + n_ox*molmass_fun('O'))';

% molm = molmass_fun(Cname);
phs_modes = zeros(length(T2d(:)),length(phs_name));
rhos      = zeros(length(T2d(:)),1);
rhof      = zeros(length(T2d(:)),1);
cwt_solid = zeros(length(Cname),length(T2d(:)));
cwt_fluid = zeros(length(Cname),length(T2d(:)));
Nsfu     = zeros(length(Cname),length(phs_name),length(T2d(:)));
for iPT = 1:length(T2d(:))
    [pc_id,phi,Cwt,Npc,rho,mu,p_out,phiw,mu2] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'CORK','S14');
    asm_id(iPT,1:length(phi)) = pc_id;
    phs_modes(iPT,pc_id)      = phi;
    if ~isempty(fluid)
        fluid_id                  =    (strcmp(phs_name(pc_id),fluid) | strcmp(phs_name(pc_id),'Melt(H18)'));        
    else
        fluid_id = [];
    end
    melt_id                   =    strcmp(phs_name(pc_id),'Melt(H18)');
    solid_id                  =   ~(strcmp(phs_name(pc_id),fluid) | strcmp(phs_name(pc_id),'Melt(H18)'));
    rhos(iPT)                 = rho(solid_id)'*phi(solid_id)/sum(phi(solid_id));
    if ~isempty(fluid_id) && sum(fluid_id)~=0
        rhof(iPT)                 = rho(fluid_id)'*phi(fluid_id)/sum(phi(fluid_id));
        cwt_fluid(:,iPT)          = Cwt(:,fluid_id)*phiw(fluid_id)/sum(Cwt(:,fluid_id)*phiw(fluid_id));
    else
        rhof(iPT) = nan;
        cwt_fluid(:,iPT) = nan;
    end
    cwt_solid(:,iPT)          = Cwt(:,solid_id)*phiw(solid_id)/sum(Cwt(:,solid_id)*phiw(solid_id));    
    Nsfu(:,pc_id,iPT)        = Npc;
end
Xh                            = cwt_solid(strcmp(Cname,'H'),:)';
[stb_phs_modes,stb_phs_names,stb_solid_modes,stb_solid_names] = PhaseModes2StablePhases(phs_modes,phs_name,fluid);
if ndim == 1
    switch ps_type
        case 'X'
            figure(1),area(squeeze(X2d),stb_solid_modes,'FaceColor','flat'),xlabel('X'),legend(stb_solid_names),axis tight            
        case 'P'
            figure(1),area(squeeze(P2d)/1e9,stb_solid_modes,'FaceColor','flat'),xlabel('P(GPa)'),legend(stb_solid_names),axis tight
            figure(2),plot(squeeze(P2d)/1e9,rhos),xlabel('P(GPa)'),title('\rho_s (kg/m^3)')
            figure(3),plot(squeeze(P2d)/1e9,Xh),xlabel('P(GPa)'),title('X_H (wt)')
            figure(4),plot(squeeze(P2d)/1e9,(1-Xh).*rhos),xlabel('P(GPa)'),title('rho_n (kg/m^3)')
        case 'T'
            figure(1),area(squeeze(T2d)-273.15,stb_solid_modes,'FaceColor','flat'),xlabel('T(\circC)'),legend(stb_solid_names),axis tight
    end
elseif ndim == 2
    switch ps_type
        case 'T-X'
            stb_phs_modes   = reshape(stb_phs_modes,length(X),length(T),length(stb_phs_names));
            stb_solid_modes = reshape(stb_solid_modes,length(X),length(T),length(stb_solid_names));            
            figure(1),
            tl_psection(X,T-273.15,Cname,asm_id,phs_name,0,[0 0],10);
            figure(2),
            nrow          = fix(sqrt(length(stb_solid_names)));% Make number of rows for subplot
            ncol          = ceil(length(stb_solid_names)/nrow);% Make number of columns for subplot
            for ipl = 1:length(stb_solid_names)
                subplot(nrow,ncol,ipl),contourf(X2d,T2d-273.15,stb_solid_modes(:,:,ipl)*100),colorbar
                title(stb_solid_names(ipl)),xlabel('T(\circC)'),xlabel('P(GPa)')
            end
        case 'P-T'
            stb_phs_modes   = reshape(stb_phs_modes,length(T),length(P),length(stb_phs_names));
            stb_solid_modes = reshape(stb_solid_modes,length(T),length(P),length(stb_solid_names));
            cwt_solid       = reshape(cwt_solid,length(Cname),length(T),length(P));
            cwt_fluid       = reshape(cwt_fluid,length(Cname),length(T),length(P));
            figure(1),
            tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name,0,[0 0],10);
            % figure(2),
            % nrow          = fix(sqrt(length(stb_solid_names)));% Make number of rows for subplot
            % ncol          = ceil(length(stb_solid_names)/nrow);% Make number of columns for subplot
            % for ipl = 1:length(stb_solid_names)
            %     subplot(nrow,ncol,ipl),contourf(T2d-273.15,P2d/1e9,stb_solid_modes(:,:,ipl)*100),colorbar
            %     title(stb_solid_names(ipl)),xlabel('T(\circC)'),xlabel('P(GPa)')
            % end
            % figure(3), 
            % rhos = reshape(rhos,length(T),length(P));
            % contourf(T2d-273.15,P2d/1e9,rhos),colorbar
            % title('\rho_s (kg/m^3)'),xlabel('T(\circC)'),xlabel('P(GPa)')      
            % figure(4), 
            % Xh = reshape(Xh,length(T),length(P));
            % contourf(T2d-273.15,P2d/1e9,Xh),colorbar
            % title('X_H (wt)'),xlabel('T(\circC)'),xlabel('P(GPa)')      
            % figure(5), 
            % rhof = reshape(rhof,length(T),length(P));
            % contourf(T2d-273.15,P2d/1e9,rhof),colorbar
            % title('\rho_f (kg/m^3)'),xlabel('T(\circC)'),xlabel('P(GPa)')      
    end
end