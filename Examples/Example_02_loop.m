clear, addpath ../ ../Utilities/ ../Solutions/ ../EOS
Temp   = linspace(300,1000,30) + 273.15;
Pres   = linspace(0.1,4,31)*1e9;
phase  = {'q,tc-ds633','coe,tc-ds633'};
for iT = 1:length(Temp)
    for iP = 1:length(Pres)
        g(iT,iP,:) = tl_gibbs_energy(Temp(iT),Pres(iP),phase);
    end
end