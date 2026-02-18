clear,close all,colormap(jet),addpath ../ ../EOS
% Input
T          = linspace(300,800,200)+ 273.15;                % Temperature (K)
P          = linspace(0.1e9,1e9,200);                      % Pressure (Pa)
% phase      = {'ky,tc-ds55','sill,tc-ds55','and,tc-ds55'};  % Phase name
phase      = {'ky,tc-ds633','sill,tc-ds633','and,tc-ds633'};  % Phase name
% Numerics
[T2d,P2d]  = ndgrid(T,P);                                  % Create P-T grid
% Gibbs energy function
g          = tl_gibbs_energy(T2d(:),P2d(:),phase);         % G (J/mol)
% Find minimum of Gibbs energy and stable phase index (ind)
[gmin,ind] = min(g); 
% Visualize stable phase
pcolor(T2d-273.15,P2d/1e9,reshape(ind,length(T),length(P))),colorbar,shading flat 