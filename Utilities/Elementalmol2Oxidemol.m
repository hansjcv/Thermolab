function [Nsys,Cname_Oxide] = Elementalmol2Oxidemol(Cname,molElems)
Elems_list  = {'Si' ,'Al' , 'Cr',    'Ti'     ,'Fe'   ,'Mn',    'Mg',    'Ca',   'Na',    'K',    'H','O'  };
Oxide_list = {'SiO2','Al2O3','Cr2O3','TiO2',  'FeO'  ,'MnO', 'MgO',  'CaO',  'Na2O', 'K2O',  'H2O'};
c_ind = zeros(size(Cname));
icat = ~strcmp(Cname,'O');
for i_c = 1:length(Cname)
    if ~strcmp(Cname(i_c),'O')
        c_ind(i_c) = find(strcmp(Elems_list,Cname(i_c)));
    end
end
c_ind(c_ind==0)=[];
Cname_Oxide = [Oxide_list(c_ind)];
noxy   = [2      3      3        2        1       1         1        1       1        1       1       ];
ncat   = [1      2      2        1        1       1         1        1       2        2       2       ];
% molmOx = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysElems = molElems(icat);
Nsys   = NsysElems./ncat(c_ind);
% Nsys   = [Nsys Nsys*(noxy(c_ind)./ncat(c_ind))']; 
Nsys   = Nsys/sum(Nsys);
