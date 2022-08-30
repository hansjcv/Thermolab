function [g,Npc,pc_id,p,z,g_id] = tl_gibbs_energy(T,P,phase,td,p,g0,v0,rho_w,eps_di)
%% tl_gibbs_energy Thermolab Gibbs energy.
%   g = tl_gibbs_energy(T,P,phase) gives the Gibbs energy of an endmember
%   at T in Kelvin and P in Pascal. Gibbs energies for multiple P-T points
%   maybe obtained by using a vector of T and a corresponding vector of P.
%   Phases must be cell array of strings holding the phase abbreviation.
%
%   [g,Npc] = tl_gibbs_energy(T,P,phase) gives additionally Npc the molar composition
%   of the phases.
%   
%   [g,Npc] = tl_gibbs_energy(T,P,phase,td) gives the Gibbs energy of a solution
%   where td is a structure array holding solution model data. Use
%   init_thermo to generate the td structure.
% 
%   [g,Npc,pc_id] = tl_gibbs_energy(T,P,phase,td) gives a unique identifier
%   for each separate phase. Phases that are discretized compounds of a
%   solutions get the same identifier
%
%   [g,Npc,pc_id,p] = tl_gibbs_energy(T,P,phase,td) outputs the proportions
%   of each solution in a cell array p.
%
%   [g,Npc,pc_id] = tl_gibbs_energy(T,P,phase,td,p) uses the user defined
%   grid of proportions of each solution in a cell array p to compute g.
%
%   Written by J. C. Vrijmoed (2022)
if ~exist('td','var')    
    td    = init_thermo(phase);% Load thermodynamic data
end
% Water properties
if ~exist('rho_w','var')
    [rho_w,eps_di] = water_props(T,P,phase);
end
if ~exist('g0','var')
    for i_sol = 1:length(phase)
        [g0(i_sol),v0(i_sol)] = tl_g0(T,P,td(i_sol),rho_w(:,i_sol),eps_di(:,i_sol));
    end
end
% Proportions and site fractions
if ~exist('p','var')
    p = props_generate(td);
end
for i_sol = 1:length(phase)
    mtpl  = td(i_sol).mtpl; zt = td(i_sol).zt; alp = td(i_sol).alp; w = td(i_sol).w; z_tol = td(i_sol).z_tol;
    z     = p{i_sol}*zt;      z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;mod_id = td(i_sol).mod_id;
    chg   = td(i_sol).chg;n_em = td(i_sol).n_em;
    if mod_id == 4        
        [mtpl,z,zt] = temkin_H18(phase(i_sol),td(i_sol),p{i_sol},z);        
    end  
    % Gibbs energy of mixing
    g_mech = p{i_sol}*g0{i_sol}';                                                                 % Mechanical
    Sconf =   8.3144*sum(mtpl.*(z.*log(z+double(z==0))- p{i_sol}*(zt.*log(zt+double(zt==0)))),2); % Configurational Entropy
    for iPT = 1:length(T)        
        g_id     = T(iPT)*Sconf;                           % Ideal mixing
        g_nid    = gnid(T(iPT),P(iPT),p{i_sol},mod_id,v0{i_sol}(iPT,:),alp,w,rho_w(iPT,i_sol),eps_di(iPT,i_sol),chg);   % Non-ideal
        g{i_sol}(:,iPT) = g_mech(:,iPT) + g_id + g_nid; % Gibbs energy in Joule/mol
    end    
    Npc{i_sol}   = p{i_sol}*n_em;
    pc_id{i_sol} = ones(1,size(p{i_sol},1))*i_sol;
end
g     = cell2mat(g');
Npc   = cell2mat(Npc')';
pc_id = cell2mat(pc_id);
end