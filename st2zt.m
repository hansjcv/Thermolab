function [zt,zmax] = st2zt(st,site_id)
zmax = zeros(size(st));
for i_site = 1:max(site_id)
    zmax(:,site_id==i_site) = repmat(sum(st(:,site_id==i_site),2),1,sum(site_id==i_site));
end
zt = st./zmax;   zt(isnan(zt)) = 0;
end