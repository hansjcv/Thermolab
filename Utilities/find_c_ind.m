function c_ind = find_c_ind(Cname,elements)
for i_c = 1:length(elements)
        c_curr = find(strcmp(Cname,elements(i_c)));
        if ~isempty(c_curr)
            c_ind(i_c) = c_curr;
        else
            c_ind(i_c) = -1;
        end
end 
end