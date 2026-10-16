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