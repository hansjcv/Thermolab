function [mtpl,z,zt] = temkin_H18(sol_id,td,p,z)
zt     = td.zt;
if strcmp(sol_id,'Melt(H18)') 
    if sum(strcmp(td.p_name,'ctL,tc-ds633'))>0
        yct   = p(:,strcmp(td.p_name,'ctL,tc-ds633'));
    else
        yct = zeros(size(p,1),1);
    end
    eval(td.z_fac)
    z = z.*z_fac;
    eval(td.z_name)
    eval(td.mtpl)
    zt    = ones(size(zt));
end
end