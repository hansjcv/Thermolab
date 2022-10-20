function [asm_legend,asm_name,asm_id] = tl_psection(X,Y,comp_names,asm_id,phs_names,tri_plot,field_size,font_size)
if nargin < 6
    tri_plot = 0;
    font_size = 10;
    min_field_size = 0;
    max_field_size = 0;
elseif nargin < 7
    font_size = 10;
    min_field_size = 0;
    max_field_size = 0;
elseif nargin < 8
    font_size = 10;
    min_field_size = field_size(1);
    max_field_size = field_size(2);
else
    min_field_size = field_size(1);
    max_field_size = field_size(2);
end
% find all unique assemblages
asm_id_raw  = sort(asm_id,2);
asm_ids     = unique(sort(asm_id_raw,2),'rows');
asm_id_arr = zeros(size(asm_id,1),1);
for i = 1:size(asm_ids)
    asm_id_arr(sum(abs(asm_id_raw-asm_ids(i,:)),2)==0) = i;
    asm_name{i} = phs_names(asm_ids(i,asm_ids(i,:)>0));
end
if tri_plot == 0
    asm_id = reshape(asm_id_arr,length(X),length(Y));
    % Make a phase diagram section with fields coloured by variance
    clf,
    field_var = zeros(size(asm_id));
    for k = 1:length(asm_name)
        field_var(asm_id==k) = length(comp_names) - length(asm_name{k});
    end
    field_var(field_var<0) = nan;
    n_var = max(field_var(:))-min(field_var(:));
    my_color = [1 1 1];
    map = 0.8:-(0.8-0.2)/n_var:0.2;
    mix = 0.5;
    map = repmat(map',1,3)*(1-mix) + repmat(my_color,size(map,2),1)*mix;
    imagesc(X,Y,field_var'),colormap(map),hold on,shading flat,set(gca,'YDir','normal')
    x = X;
    y = Y;
    z = asm_id;
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    % Brute force way to draw lines around phase fields
    hold on
    for i = 1:length(x)
        for j = 1:length(y)
            % bottom
            if j > 1
                x_line = [x(i)-dx/2 x(i)+dx/2];
                y_line = [y(j)-dy/2 y(j)-dy/2];
                if z(i,j-1)~=z(i,j)
                    plot(x_line',y_line','k-','LineWidth',1)
                end
            end
            % right
            if i < length(x)
                x_line = [x(i)+dx/2 x(i)+dx/2];
                y_line = [y(j)-dy/2 y(j)+dy/2];
                if z(i,j)~=z(i+1,j)
                    plot(x_line',y_line','k-','LineWidth',1)
                end
            end
            % top
            if j < length(y)
                x_line = [x(i)-dx/2 x(i)+dx/2];
                y_line = [y(j)+dy/2 y(j)+dy/2];
                if z(i,j)~=z(i,j+1)
                    plot(x_line',y_line','k-','LineWidth',1)
                end
            end
            % left
            if i > 1
                x_line = [x(i)-dx/2 x(i)-dx/2];
                y_line = [y(j)-dy/2 y(j)+dy/2];
                if z(i-1,j)~=z(i,j)
                    plot(x_line',y_line','k-','LineWidth',1)
                end
            end
        end
    end
    [X2d,Y2d] = ndgrid(X,Y);
else
    asm_id = reshape(asm_id_arr,size(X,1),size(Y,2));    
    % Make a phase diagram section with fields coloured by variance
    clf,
    field_var = nan(size(asm_id));
    for k = 1:length(asm_name)
        field_var(asm_id==k) = length(comp_names) - length(asm_name{k});
    end
    field_var(field_var<0) = nan;
    n_var = max(field_var(:))-min(field_var(:));
    my_color = [1 1 1];
    map = 0.8:-(0.8-0.2)/n_var:0.2;
    mix = 0.5;
    map = repmat(map',1,3)*(1-mix) + repmat(my_color,size(map,2),1)*mix;    
    for i_field = 1:max(asm_id(:))
        c{i_field} = contour(X,Y,asm_id==i_field,[0.5,0.5],'k','Visible','off');
    end
    pcolor(X,Y,field_var),colormap(map),hold on,shading flat,set(gca,'YDir','normal'),axis equal
    plot([0,1/sind(60);1/sind(60),1/sind(60)/2;1/sind(60)/2,0],[0,0;0,1;1,0],'k'),axis([0 1/sind(60) 0 1])    
    for i_field = 1:max(asm_id(:))
        plot(c{i_field}(1,2:c{i_field}(2,1)),c{i_field}(2,2:c{i_field}(2,1)),'k');
        hold on
    end
    X2d = X;
    Y2d = Y;
end
% Add field labels
cnt = 0;
asm_legend = [];
for j = 1:length(asm_name)
    if sum(asm_id(:)==j)>max_field_size
        t_h = text(mean(X2d(asm_id==j)),mean(Y2d(asm_id==j)),asm_name{j},'Color','k','BackgroundColor','w','HorizontalAlignment','center','EdgeColor','k','FontSize',font_size);
    elseif sum(asm_id(:)==j)>min_field_size && sum(asm_id(:)==j)<max_field_size
        cnt = cnt + 1;
        asm_legend{cnt} = asm_name{j};
        text(mean(X2d(asm_id==j)),mean(Y2d(asm_id==j)),num2str(cnt),'Color','k','BackgroundColor','w','HorizontalAlignment','center','EdgeColor','k','FontSize',font_size);
    end
end