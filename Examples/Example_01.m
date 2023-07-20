clear, addpath ../
T     = 300 + 273.15;
P     = 1e9;
phase = {'q,tc-ds633'};
g     = tl_gibbs_energy(T,P,phase);