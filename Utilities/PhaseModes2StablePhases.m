function [stb_phs_modes,stb_phs_names,stb_solid_modes,stb_solid_names] = PhaseModes2StablePhases(phs_modes,phs_name,fluid)
stb_phs_modes   = phs_modes(:,sum(phs_modes,1)>0);
stb_phs_names   = phs_name(sum(phs_modes,1)>0);
stb_solid_modes = stb_phs_modes(:,~strcmp(stb_phs_names,fluid));
stb_solid_names = stb_phs_names(~strcmp(stb_phs_names,fluid));
stb_solid_modes = stb_solid_modes./repmat(sum(stb_solid_modes,2),1,length(stb_solid_names)+1e-20);
% stb_phs_modes   = reshape(stb_phs_modes,length(T),length(P),length(stb_phs_names));
% stb_solid_modes = reshape(stb_solid_modes,length(T),length(P),length(stb_solid_names));