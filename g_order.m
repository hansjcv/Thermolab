function [g,p_all] = g_order(T,P,phs_name,td,c1,z_od,g0,V0)
z_tol = 1e-10;
[c1_arr,z_od_arr] = ndgrid(c1,z_od);
p{1} = (td.p_from_c_cons*[ones(1,length(c1_arr(:)));c1_arr(:)';z_od_arr(:)'])';
z     = p{1}*td.zt;
z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;
badz = min(z')<0|max(z')>1;
p{1}(badz,:) = [];
p{1}(p{1}>1-z_tol&p{1}<1+z_tol) = 1;p{1}(p{1}>0-z_tol&p{1}<0+z_tol) = 0;
[T2d,P2d] = ndgrid(T,P);
if ~exist('g0','var')
    [g0,V0]  = tl_g0(T2d(:),P2d(:),td); % get endmember gibbs energies of the solution model
end
g_full = tl_gibbs_energy(T2d(:),P2d(:),phs_name,td,p,g0,V0);
g_all = nan(length(c1),length(z_od));
for iPT = 1:length(T2d(:))
    g_all(~badz) = g_full(:,iPT);
    [gmin,id] = min(g_all,[],2);
    g(:,iPT) = gmin;
    p{1} = (td.p_from_c_cons*[ones(1,length(c1));c1;z_od(id)])';
    p{1}(p{1}>1-z_tol&p{1}<1+z_tol) = 1;p{1}(p{1}>0-z_tol&p{1}<0+z_tol) = 0;
    p_all{iPT} = p;
end
