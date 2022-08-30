function phases = get_phases(db_name,solfile,Cname)
% clear,clf
addpath ../
load('tl_dataset','elements','db_names','phs_names','nphs')
% db_name = 'tc-ds633';
% solfile = 'solution_models_H18';
% Cname   = {'Si','Mg','Fe','O'};
% Cname   ={'Si' ,'Al'    'Fe'   ,'Mg',    'H',   'K',     'O' }
phases = phs_names(strcmp(db_names,db_name));
for i_c = 1:length(Cname)
    c_ind(i_c) = find(strcmp(elements,Cname(i_c)));
end
c_excl        = ones(1,length(elements));
c_excl(c_ind) = 0;
phases = phases(sum(abs(nphs(strcmp(db_names,db_name),c_excl==1)),2)==0&sum(abs(nphs(strcmp(db_names,db_name),c_excl==0)),2)>0);
if exist('sheetnames')==2
    sol_mod_list = sheetnames([solfile '.xlsx']);
else
    [status,sol_mod_list] = xlsfinfo([solfile '.xlsx']);
end
% phases = [];
for i_sol = 1:length(sol_mod_list)
    td = init_thermo(sol_mod_list(i_sol),elements,solfile);
    if sum(sum(abs(td.n_em(:,c_excl==1)),2)==0&sum(abs(td.n_em(:,c_excl==0)),2)>0)>0
        phases = [phases sol_mod_list(i_sol)];
    end
    for ip = 1:length(td.p_name)
        if sum(strcmp(phases,td.p_name(ip)))>0
            phases(strcmp(phases,td.p_name(ip))) = [];
        end
    end
end