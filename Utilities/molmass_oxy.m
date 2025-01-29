function molm = molmass_oxy(Cname)
n_ox = zeros(size(Cname));
for i = 1:length(Cname)
    ox_string = Cname{i};
    digit_id  = isstrprop(ox_string,'digit');
    ox_letter = ox_string(~digit_id);
    upper_id  = find(isstrprop(ox_letter,'upper'));
    if numel(upper_id)>1
        cat_name{i} = ox_letter(upper_id(1):upper_id(2)-1);
    else
        cat_name{i} = ox_letter(upper_id(1):end);
        n_cat(i) = 1;
        continue
    end
    nele_id     = find(digit_id);   
    if sum(digit_id)==0 % in case no numbers are found
        n_cat(i)    = 1;
        n_ox(i)     = 1;
    elseif numel(nele_id)==2
        n_cat(i)    = eval(ox_string(nele_id(1)));
        n_ox(i)     = eval(ox_string(nele_id(2)));
    elseif numel(nele_id)==1 && strfind(Cname{i},'O')~=numel(ox_string) % now there are more than one O's
        n_cat(i)    = 1;
        n_ox(i)     = eval(ox_string(nele_id(1)));
    elseif numel(nele_id)==1 && strfind(Cname{i},'O')==numel(ox_string) % now there are more than one O's
        n_cat(i)    = eval(ox_string(nele_id(1)));
        n_ox(i)     = 1;
    end
end
molm = (n_cat.*molmass_fun(cat_name)' + n_ox*molmass_fun('O'))';