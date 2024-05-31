function disp_reactions(phases,v)
n_react = size(v,2);
for i_react = 1:n_react
    reactants             = phases(v(:,i_react)<0);
    v_reactants           = v(v(:,i_react)<0,i_react);
    products              = phases(v(:,i_react)>0);
    v_products            = v(v(:,i_react)>0,i_react);
    reaction_string       = [];
    for ircnt = 1:length(reactants)
        if ircnt < length(reactants)
            reaction_string = [reaction_string num2str(-v_reactants(ircnt)) ' ' reactants{ircnt} ' + '];
        else
            reaction_string = [reaction_string num2str(-v_reactants(ircnt)) ' ' reactants{ircnt} ' = '];
        end
    end
    for ipcnt = 1:length(products)
        if ipcnt < length(products)
            reaction_string = [reaction_string num2str(v_products(ipcnt)) ' ' products{ipcnt} ' + '];
        else
            reaction_string = [reaction_string num2str(v_products(ipcnt)) ' ' products{ipcnt}];
        end
    end
    disp(reaction_string)
end