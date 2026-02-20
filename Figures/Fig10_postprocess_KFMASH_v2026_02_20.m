clear,addpath ../ ../Utilities/
runname = 'KFMASH_2026_02_20_80x81';
load(['linprog_run_' runname]);                                          % load linprog run data
molm = molmass_fun(Cname);
solv_tol = 1;
fluid = 'H2O,tc-ds55';
[rhos,rhof,cwt_solid,cwt_fluid,phs_modes,Nsfu,asm_id] = postprocess_results(runname,fluid,2,'P-T',solv_tol);