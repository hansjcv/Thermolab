clear,clf
ele                = {'Si','Al','O'};
unknown_min(1,1,:) = [1     0     2];
unknown_min(2,1,:) = [1     3     5];
% List of minerals
min_name = {'Quartz','Alum'};
comp     = {'Si','Al','O'};
sfu      = [1    0  2
            1    2  5];
for i = 1:length(comp)
    c_ind(i) = find(strcmp(ele,comp(i)));
end
ind = zeros(2,length(ele),length(min_name));
dc = 0.001;
for k = 1:length(min_name)
    ind(1,c_ind,k) = sfu(k,:)-dc;
    ind(2,c_ind,k) = sfu(k,:)+dc;
end
C_maps_all = unknown_min;
%% Recognizing phases
phase_map = zeros(size(C_maps_all,1),size(C_maps_all,2));
for ip = 1:size(ind,3)
    find_id = ['id_ic = C_maps_all(:,:,' num2str(1) ')>=ind(1,' num2str(1) ',' num2str(ip) ') & C_maps_all(:,:,' num2str(1) ')<=ind(2,' num2str(1)  ',' num2str(ip) ')'];
    for ic = 2:length(ele)
        find_id = [find_id '& C_maps_all(:,:,' num2str(ic) ')>=ind(1,' num2str(ic)  ',' num2str(ip) ') & C_maps_all(:,:,' num2str(ic) ')<=ind(2,' num2str(ic) ',' num2str(ip) ')'];
    end
    if ip == 1
        eval([find_id ';']);
    else
        eval([find_id '& phase_map==0;']);
    end
    phase_map(id_ic)=ip;    
end
min_name = ['unknown',min_name]
phases = min_name(phase_map+1)