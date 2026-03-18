function [p,z_indep] = props_generate(td)
for i_sol = 1:length(td)
    z_indep{i_sol}  = auto_xgrid(td(i_sol).z_lim,td(i_sol).nc,td(i_sol).subdtype);
    if isempty(td(i_sol).site_var)
        p{i_sol} = 1;
        z_indep{i_sol} = 1;
    else                
        if td(i_sol).mod_id == 4
            [p{i_sol},z_indep{i_sol}] = zvar2zap(z_indep{i_sol},td(i_sol).phase_name,td(i_sol).site_id,td(i_sol).site_var,td(i_sol).p_from_z_cons,td(i_sol).st,td(i_sol).z_tol);
        else
            [p{i_sol},z_indep{i_sol}] = zvar2zap(z_indep{i_sol},td(i_sol).phase_name,td(i_sol).site_id,td(i_sol).site_var,td(i_sol).p_from_z_cons,td(i_sol).zt,td(i_sol).z_tol);
        end
    end
end
end
function [x_arr,x] = auto_xgrid(x_in,nc,subd_type)
if ~isempty(x_in)    
    if nargin < 3
        subd_type = zeros(1,size(x_in,2));
    end
    nvar = size(x_in,2);
    if isscalar(nc)
        nc = ones(1,nvar)*nc;
    end
    for i = 1:size(x_in,2)
        if subd_type(i) ==1
            x_in(:,i) = log10(x_in(:,i)+double(x_in(:,i)==0));
        end
    end
    n_grid = '[';
    for i_var = 1:nvar
        if i_var<nvar
            n_grid = [n_grid 'x' num2str(i_var) ','];
        else
            n_grid = [n_grid 'x' num2str(i_var) '] = ndgrid('];
        end
    end
    ngrid_arr = 'x_arr = [';
    for i_var = 1:nvar
        if i_var<nvar
            if subd_type(i_var) == 0
                n_grid = [n_grid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')),'];
            elseif subd_type(i_var) == 1
                if x_in(1,i_var)==0
                    x_in(1,i_var) = -16;
                    n_grid = [n_grid '[0,logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) '))],'];
                else
                    n_grid = [n_grid 'logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')),'];
                end
            elseif subd_type(i_var) > 1
                n_grid = [n_grid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')).^' num2str(subd_type(i_var)) ','];
            end
            ngrid_arr = [ngrid_arr 'x' num2str(i_var) '(:),'];
        else
            if subd_type(i_var) == 0
                n_grid = [n_grid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')));'];
            elseif subd_type(i_var) == 1
                if x_in(1,i_var)==0
                    x_in(1,i_var) = -16;
                    n_grid = [n_grid '[0,logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) '))]);'];
                else
                    n_grid = [n_grid 'logspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')));'];
                end
            elseif subd_type(i_var) > 1 
                n_grid = [n_grid 'linspace(x_in(1,' num2str(i_var) '),x_in(2,' num2str(i_var) '),nc(' num2str(i_var) ')).^' num2str(subd_type(i_var)) ');'];
            end
            ngrid_arr = [ngrid_arr 'x' num2str(i_var) '(:)];'];
        end
    end
    eval(n_grid);
    eval(ngrid_arr);
    for i_var = 1:nvar
        x{i_var} = eval(['x' num2str(i_var) ';']);
    end
else
    x_arr = 1;
    x{1} = 1;
end
end
function [p,z_indep] = zvar2zap(z_indep,sol_id,site_id,site_var,p_from_z_cons,zt,z_tol)
indep_site_id = site_id(site_var);
if ~strcmp(sol_id,'Melt(G25)')
    for i_site = 1:max(site_id)
        z_indep(sum(z_indep(:,indep_site_id==i_site),2)>1,:)=[];
    end
end
% Site fractions and proportions
p    = (p_from_z_cons*[ones(1,size(z_indep,1)); z_indep'])';
z    = p*zt;
z(z<1+z_tol & z>1-z_tol) = 1;z(z<  z_tol & z> -z_tol) = 1e-20;
if ~strcmp(sol_id,'Melt(G25)')
    badz = min(z')<0|max(z')>1;
else
    badz = min(z')<0;
end
p(badz,:)       = [];
z_indep(badz,:) = [];
end