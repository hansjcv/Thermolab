function [Nsys,Cname] = Oxidemol2Elementalmol(Cname_Oxide,molOx)
Elems_list  = {'Si' ,'Al' , 'Cr',    'Ti'     ,'Fe'   ,'Mn',    'Mg',    'Ca',   'Na',    'K', 'C',   'H',  'O'  };
Oxide_list = {'SiO2','Al2O3','Cr2O3','TiO2',  'FeO'  ,'MnO',    'MgO',  'CaO',  'Na2O', 'K2O', 'CO2', 'H2O' };
c_ind = zeros(size(Cname_Oxide));
O_ind = zeros(size(Cname_Oxide));
for i_c = 1:length(Cname_Oxide)
    if ~strcmp('O',Cname_Oxide(i_c))
        c_ind(i_c) = find(strcmp(Oxide_list,Cname_Oxide(i_c)));
    else
        O_ind(i_c) = 1;%find(strcmp(Elems_list,'O'));
    end
end
Cname = [Elems_list(c_ind(c_ind~=0)),'O'];
noxy   = [2      3      3        2        1       1         1        1       1        1     2    1       ];
ncat   = [1      2      2        1        1       1         1        1       2        2     1    2       ];
% molmOx = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysOx = molOx(~O_ind);
Nsys   = NsysOx.*ncat(c_ind(c_ind~=0));
Nsys   = [Nsys Nsys*(noxy(c_ind(c_ind~=0))./ncat(c_ind(c_ind~=0)))']; 
if sum(strcmp(Cname_Oxide,'O'))==1
    Nsys(end) = Nsys(end) + molOx(O_ind==1);
end
Nsys   = Nsys/sum(Nsys);
