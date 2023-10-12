function compound_wt = molmass_fun(compounds,compound_atoms,nmol)
%% Function to acquire weights of compounds (or phases)
% Updated 31.08.2023, NOTE: elements with zero mass have no standard atomic
% mass as they are never stable isotopes
% Usage:
% compound_wt = atomwt_fun(compounds,compound_atoms,nmol)
%% Example input
% compounds      = {'fo' 'en' 'per' 'q'};
% compound_atoms = {'Mg' 'Si' 'O'};
% nmol        = [ 2     1   4
%                 2     2   6 
%                 1     0   1
%                 0     1   2];
if nargin < 2
    nmol = eye(length(compounds),length(compounds));
    compound_atoms = compounds;
end
load atomweight atomwt5 symbol
atomwt5 = [0; atomwt5];
atomwt5(isnan(atomwt5)) = 0;
symbol{87} = 'Fr';
symbol  =['e'; symbol];
wt = zeros(size(nmol));
for i = 1:length(compounds)    
    for j = 1:length(compound_atoms)        
        wt(i,j) = nmol(i,j)*atomwt5(strcmpi(symbol,compound_atoms(j)));
    end
end
compound_wt = sum(wt,2)/1e3;
end