function [mu,a,RTlngam,z] = tl_chemical_potential(T,P,td,p,g0,v0)
R = 8.3144;P = P/1e8;
st = td.st;zt = td.zt;mtpl = td.mtpl;alp = td.alp;w = td.w;mod_id = td.mod_id;z_tol = td.z_tol;
if mod_id == 0 || mod_id == 1 || mod_id == 4 || mod_id == 5
    if mod_id == 4
        zt = st;
    end   
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
elseif mod_id == 3
    V = cell2mat(v0);
    % Mixing parameters Aranovich & Newton 1999
    w_data   = [12893, -6.501,1.0112];
    A = w_data(1);
    B = w_data(2);
    C = w_data(3);
    W_AN99 = (A+B*T).*(1-exp(-20*P)) + C*T.*P;
    RTlngam(:,1) = p(:,2).^2*W_AN99.*(V(1).*V(2).^2./((V(1)+V(2)).*(p(:,1)*V(1)+p(:,2)*V(2)).^2));
    RTlngam(:,2) = p(:,1).^2*W_AN99.*(V(2).*V(1).^2./((V(1)+V(2)).*(p(:,1)*V(1)+p(:,2)*V(2)).^2));
    a = p; 
    a(a==0) = z_tol;
end
mu0 = cell2mat(g0);
if mod_id == 5
    mu0(1:end-1) = mu0(1:end-1) + R*T*log(1000/18.0150);
end
mu = mu0 + R*T*log(a) + RTlngam;

