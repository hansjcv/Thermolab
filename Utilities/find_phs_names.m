function phs_db_name = find_phs_names(db_name,element)
load tl_dataset;
nphs_db_name = nphs(contains(phs_names,db_name),:); 
phs_db_name = phs_names(contains(phs_names,db_name));
phs_db_name = phs_db_name(nphs_db_name(:,strcmp(elements,element))>0)';