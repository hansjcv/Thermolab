function [Nsys,Cname] = Oxideweight2Elementalmol(Cname_Oxide,wtOx)
Elems_list = {'Si'  ,'Al'   , 'Cr',    'Ti','Fe'  , 'Fe'     ,'Mn',    'Mg',    'Ca',    'Na',    'K',    'H','O'  };
Oxide_list = {'SiO2','Al2O3','Cr2O3','TiO2','FeO' , 'Fe2O3'  ,'MnO', 'MgO' ,  'CaO',  'Na2O', 'K2O',  'H2O'};
noxy       = [2     ,3      ,3      , 2    , 1    , 3        , 1      ,   1  ,     1,       1,     1,     1       ];
ncat       = [1     ,2      ,2      , 1    , 1    , 2        , 1      ,   1  ,     1,       2,     2,     2       ];
molmOx     = [60.084,101.961,151.9904, 79.8658,71.844, 159.6882 ,70.93744,40.304, 56.077, 61.97894, 94.196,  18.01528];
c_ind = zeros(size(Cname_Oxide));
for i_c = 1:length(Cname_Oxide)
    c_ind(i_c) = find(strcmp(Oxide_list,Cname_Oxide(i_c)));
end
Cname_all = [Elems_list(c_ind),'O'];
NsysOx    = wtOx./molmOx(c_ind);
Nsys_all  = NsysOx.*ncat(c_ind);
Nsys_all  = [Nsys_all Nsys_all*(noxy(c_ind)./ncat(c_ind))'];
Cname      = unique(Cname_all,'stable');
Nsys      = zeros(size(Cname));
for i_c = 1:length(Cname)
    Nsys(i_c) = sum(Nsys_all(strcmp(Cname_all,Cname(i_c))));
end
Nsys   = Nsys/sum(Nsys);