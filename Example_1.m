clear,clf
T         = linspace(200,1000,10)+273.15;                         % Temperature in K
P         = linspace(0.1,4,10)*1e9;                               % Pressure in Pa
phase     = {'an,tc-ds55','ky,tc-ds55','gr,tc-ds55','q,tc-ds55'}; % Names of phase (HP98 dataset = tc-ds55, HP11 dataset = tc-ds633)
[T2d,P2d] = ndgrid(T,P);                                          % Make T and P matrix
[g,N]     = tl_gibbs_energy(T2d(:),P2d(:),phase);                 % Get Gibbs energies (see Thermolab paper, when published)
v         = null(N);                                              % Reaction stoichiometry
dg        = v'*g;                                                 % Delta G of reaction
% Visualize equilibrium
contourf(T2d-273.15,P2d/1e9,reshape(dg,length(T),length(P)),[0,0]);