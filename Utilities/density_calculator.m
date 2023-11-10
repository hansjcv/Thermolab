clear,clf, addpath ../
T = linspace(400,1000,100) + 273.15;
P = linspace(1e9,4e9,101);
phs_name = {'atg,tc-ds55','br,tc-ds55','fo,tc-ds55','H2O,tc-ds55'};
load tl_dataset
molm       = molmass_fun(elements);
[Temp,Pres] = ndgrid(T,P);
rho = zeros(length(Temp(:)),length(phs_name));
for ip = 1:length(phs_name)
    id    = find(strcmp(phs_names,phs_name(ip)));
    td  = td_data{id};
    [VdP,V] = intVdP(Temp(:),Pres(:),td,eos(id));
    molweight = nphs(id,:)*molm;
    rho(:,ip) = molweight./V*1e5;
end