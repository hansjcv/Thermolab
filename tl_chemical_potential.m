function [mu,a,RTlngam,z] = tl_chemical_potential(T,P,td,p,g0)
R = 8.3144;P = P/1e8;
st = td.st;zt = td.zt;mtpl = td.mtpl;alp = td.alp;w = td.w;
if td.mod_id == 4
    zt = st;
end
z_tol = 1e-10;
z = p*zt;z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;
% Compute activity normalization factors
a = zeros(size(p));
for ip = 1:size(st,1)
    a0(ip) = prod(zt(ip,st(ip,:)>0).^(zt(ip,st(ip,:)>0).*mtpl(st(ip,:)>0)));
end
N = 1./a0;
% Ideal activities
for ip = 1:size(st,1)
    a(:,ip) = N(ip)*prod(z.^repmat(zt(ip,:).*mtpl,size(z,1),1),2);
end
a(a==0) = z_tol;
% non-ideal
alp = [1 T P]*alp;
W     = w(:,:,1)   +   w(:,:,2)*T +   w(:,:,3)*P;
np      = size(p,2);
phi     = zeros(size(p));
RTlngam = zeros(size(p));
q       = zeros(size(p));
for ip = 1:np
    for i = 1:np
        phi(:,i) = p(:,i).*alp(i)./(p*alp');
        if i==ip
            q(:,i) = 1 - phi(:,i);
        else
            q(:,i) =  - phi(:,i);
        end
        for j = 1:np
            phi(:,j) = p(:,j).*alp(j)./(p*alp');
            if j==ip
                q(:,j) = 1 - phi(:,j);
            else
                q(:,j) =   - phi(:,j);
            end
            RTlngam(:,ip) = RTlngam(:,ip) - q(:,i).*q(:,j).*W(i,j).*(alp(ip)./(alp(i)+alp(j)));
        end
    end
end
mu = (cell2mat(g0) + R*T*log(a) + RTlngam);