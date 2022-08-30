function [p,z_indep] = props_generate(td)
for i_sol = 1:length(td)
    z_indep{i_sol}  = auto_xgrid(td(i_sol).z_lim,td(i_sol).dz,td(i_sol).subdtype);
    if isempty(td(i_sol).site_var)
        p{i_sol} = 1;
    else
        p{i_sol} = zvar2zap(z_indep{i_sol},td(i_sol).phase_name,td(i_sol).site_id,td(i_sol).site_var,td(i_sol).p_from_z_cons,td(i_sol).zt,td(i_sol).z_tol);
    end
end
end
function [x_arr,x] = auto_xgrid(x_in,dc,subd_type)
if ~isempty(x_in)    
    if nargin < 3
        subd_type = zeros(1,size(x_in,2));
    end
    nvar = size(x_in,2);
    if length(dc) == 1
        dc = ones(1,nvar)*dc;
    end
    for i = 1:size(x_in,2)
        nc(i) = fix((x_in(2,i)-x_in(1,i))/dc(i)+1);
        if subd_type(i) ==1
            x_in(:,i) = log10(x_in(:,i)+double(x_in(:,i)==0)*1e-16);
        end
    end
    ngrid = '[';
    for i_var = 1:nvar
        if i_var<nvar
            ngrid = [ngrid 'x' num2str(i_var) ','];
        else
            ngrid = [ngrid 'x' num2str(i_var) '] = ndgrid('];
        end
    end
    ngrid_arr = 'x_arr = [';
    for i_var = 1:nvar
        if i_var<nvar
            if subd_type(i_var) == 0
                ngrid = [ngrid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')),'];
            elseif subd_type(i_var) == 1
                ngrid = [ngrid 'logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')),'];
            end
            ngrid_arr = [ngrid_arr 'x' num2str(i_var) '(:),'];
        else
            if subd_type(i_var) == 0
                ngrid = [ngrid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')));'];
            elseif subd_type(i_var) == 1
                ngrid = [ngrid 'logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')));'];
            end
            ngrid_arr = [ngrid_arr 'x' num2str(i_var) '(:)];'];
        end
    end
    eval(ngrid);
    eval(ngrid_arr);
    for i_var = 1:nvar
        x{i_var} = eval(['x' num2str(i_var) ';']);
    end
else
    x_arr = 1;
    x{1} = 1;
end
end
function p = zvar2zap(z_indep,sol_id,site_id,site_var,p_from_z_cons,zt,z_tol)
indep_site_id = site_id(site_var);
if ~strcmp(sol_id,'Melt(H18)')
    for i_site = 1:max(site_id)
        z_indep(sum(z_indep(:,indep_site_id==i_site),2)>1,:)=[];
    end
end
% Site fractions and proportions
p    = (p_from_z_cons*[ones(1,size(z_indep,1)); z_indep'])';
z    = p*zt;
z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;
if ~strcmp(sol_id,'Melt(H18)')
    badz = min(z')<0|max(z')>1;
else
    badz = min(z')<0;
end
p(badz,:) = [];
end