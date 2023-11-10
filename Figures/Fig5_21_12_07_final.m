clear,close all,addpath ../ ../Solutions ../EOS
runname = 'example_soap';
T       = 300 + 273.15;                                                  % T in Kelvin
P       = 0.3e9;                                                         % P in Pascal
solfile = 'solution_models_soapstone';                                   % Name of solution model file
Cname   = {'Si'  ,'Al'  ,'Mg'  ,'Fe',  'Ca',  'C',   'H'   ,'O'   ,'e'}; % Components in the system
Xsys    = [0.7392;0.0451;1.0000;0.0644;0.0113;0.0232;1.3895;3.3629; 0];  % moles of each component in system
fluid   = 'Fluid-CO2-H2O';                                          % Name of fluid model
phases  = {fluid,'Chlorite','Antigorite','Clinopyroxene','Talc','Magnesite','Dolomite','Brucite','Olivine',...
          'per,tc-ds55','q,tc-ds55','ky,tc-ds55','sill,tc-ds55','and,tc-ds55'}; % Names of phases to consider
td      = init_thermo(phases,Cname,solfile);                                    % Initialize phase data
p       = props_generate(td);                                                   % Grid of compositions for mixtures
[rho_w,eps_w] = water_props(T,P,phases,'ZD05','S14');                                  % Density and Dielectric constant
[g0,v0] = tl_g0(T,P,td,rho_w,eps_w);                                     % Endmember g and v in mixture
% Compute Gibbs energy 
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phases,td,p,g0,v0,rho_w,eps_w);      % Call Gibbs energy function
% Gibbs minimization
LB     = zeros(1,length(g));                                             % stable phase amount cannot be negative.
alph   = linprog(g,[],[],Npc,Xsys,LB);                                   % The Gibbs energy minimization
% Save results
save([runname '_linprog']);