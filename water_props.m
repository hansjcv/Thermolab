function [rho_w,eps_di] = water_props(T,P,phase,rho_w_model,eps_di_model)
rho_w  = zeros(length(T),length(phase));
eps_di = zeros(length(T),length(phase));
for i_sol = 1:length(phase)
    if ~exist('rho_w_model','var')
        if contains(phase(i_sol),'supcrt') || contains(phase(i_sol),'Miron')
            rho_w_model_2  = 'JN91';
        elseif contains(phase(i_sol),'DEW')
            rho_w_model_2  = 'ZD05';
        elseif contains(phase(i_sol),'tc-ds633') || contains(phase(i_sol),'tc-ds62')            
            rho_w_model_2  = 'PS94';
        else            
            rho_w_model_2  = 'CORK';
        end
    else
        rho_w_model_2 = rho_w_model;
    end
    if ~exist('eps_di_model','var')
        if contains(phase(i_sol),'supcrt') || contains(phase(i_sol),'Miron')            
            eps_di_model_2 = 'JN91';
        elseif contains(phase(i_sol),'DEW')            
            eps_di_model_2 = 'S14';
        else%if contains(phase(i_sol),'tc-ds55') || contains(phase(i_sol),'tc-ds633')            
            eps_di_model_2 = 'S14';
        end
    else
        eps_di_model_2 = eps_di_model;
    end
    rho_w(:,i_sol)  = rho_H2O(T,P,rho_w_model_2);       % Water density
    eps_di(:,i_sol) = eps_H2O(T,P,rho_w(:,i_sol),eps_di_model_2);  % Water dielectric constant
end