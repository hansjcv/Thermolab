function [rho_w,eps_di] = water_props(T,P,phase,rho_w_model,eps_di_model)
rho_w  = zeros(length(T),length(phase));
eps_di = zeros(length(T),length(phase));
for i_sol = 1:length(phase)
    if ~exist('rho_w_model','var')
        if contains(phase(i_sol),'supcrt') || contains(phase(i_sol),'Miron')
            rho_w_model  = 'JN91';
        elseif contains(phase(i_sol),'DEW')
            rho_w_model  = 'ZD05';
        else%if contains(phase(i_sol),'tc-ds55') || contains(phase(i_sol),'tc-ds633')
            rho_w_model  = 'CORK';
        end
    end
    if ~exist('eps_di_model','var')
        if contains(phase(i_sol),'supcrt') || contains(phase(i_sol),'Miron')            
            eps_di_model = 'JN91';
        elseif contains(phase(i_sol),'DEW')            
            eps_di_model = 'S14';
        else%if contains(phase(i_sol),'tc-ds55') || contains(phase(i_sol),'tc-ds633')            
            eps_di_model = 'S14';
        end
    end
    rho_w(:,i_sol)  = rho_H2O(T,P,rho_w_model);       % Water density
    eps_di(:,i_sol) = eps_H2O(T,P,rho_w(:,i_sol),eps_di_model);  % Water dielectric constant
end