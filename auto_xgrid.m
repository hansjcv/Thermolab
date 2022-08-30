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