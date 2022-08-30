function [p_id,st,site_id,mtpl,mod_id,alp,w,Nphs,chg,td,dGex,EOS,CEOS,mcoef,Gref] = init_phase(sol_model,mineral,Cname)
load 'tl_dataset' nphs elements phs_names make_coeff td_data ceos eos dgex Gr
for ic = 1:length(Cname),c_ind(ic) = find(strcmp(elements,Cname(ic)));end
sol_mod_list = sheetnames([sol_model '.xlsx']);
if sum(strcmp(sol_mod_list,mineral))>0
    % Read data from excel tables
    [data,txt]  = xlsread([sol_model '.xlsx'],mineral,'','basic');
    mod_id      = data(1,1);
    sitenames   = txt(strcmpi(txt(:,1),'Sitenames'),2:end);
    occ_id      = find(strcmpi(txt(:,1),'Occupancy'));
    occupancy   = txt(occ_id,2:end);
    mtpl_id     = find(strcmpi(txt(:,1),'Multiplicity'));
    st          = data(occ_id+2:mtpl_id-2,1:length(sitenames));
    mtpl        = data(mtpl_id,1:length(sitenames));
    z_min           = data(strcmpi(txt(:,1),'z_min'),:);
    z_max           = data(strcmpi(txt(:,1),'z_max'),:);
    dz              = data(strcmpi(txt(:,1),'dz'),:);
    p_name          =  txt(occ_id+2:mtpl_id-2,1);
    subdtype        = data(strcmpi(txt(:,1),'subdivision'),:);
    w0_data_start   = find(strcmpi(txt(:,1),'w0'))+1;
    if ~isempty(w0_data_start)
        w0              = data(w0_data_start:w0_data_start+length(p_name)-1,1:length(p_name));
    else
        w0 = zeros(length(p_name),length(p_name));
    end
    wT_data_start   = find(strcmpi(txt(:,1),'wT'))+1;
    if ~isempty(wT_data_start)
        wT              = data(wT_data_start:wT_data_start+length(p_name)-1,1:length(p_name));
    else
        wT = zeros(size(w0));
    end
    wP_data_start   = find(strcmpi(txt(:,1),'wP'))+1;
    if ~isempty(wP_data_start)
        wP              = data(wP_data_start:wP_data_start+length(p_name)-1,1:length(p_name));
    else
        wP = zeros(size(w0));
    end
    alp_data_start  = find(strcmpi(txt(:,1),'alp'))+1;
    if ~isempty(alp_data_start)
        alp             = data(alp_data_start:alp_data_start+2,1:length(p_name));
    else
        alp = [ones(1,length(p_name));zeros(2,length(p_name))];
    end
    w(:,:,1) = w0;
    w(:,:,2) = wT;
    w(:,:,3) = wP;
    z_lim            = [z_min;z_max];
    % Extract information from data
    sites = unique(sitenames,'stable');
    site_id = zeros(1,length(sitenames'));
    for i_site = 1:length(sites)
        site_id(strcmp(sitenames,sites(i_site))) = i_site;     % index for the site (rather than a string with names)
    end       
else
    sitenames = [];
    occupancy = 'pure';
    sites = [];
    p_name = {mineral};
    alp = [
        1.000
        0.000
        0.000
        ];
    w(:,:,1)      = 0;
    w(:,:,1)      = 0;
    w(:,:,1)      = 0;
    st      = 1;
    site_id = 1;
    mtpl    = 1;
    mod_id = 0;
    z_lim      = [0;1];
    dz         = 1;
end
for k = 1:length(p_name)
    ip  = find(strcmp(phs_names,p_name(k)));                 % Find index of the endmember in the list
    Nphs(k,:) = make_coeff(ip)*nphs(ip,c_ind);
    chg(k)    = make_coeff(ip)*nphs(ip,strcmp(elements,'e'));
    p_id{k} = ip;
    td{k}   = cell2mat(td_data(ip)');
    dGex{k}   = dgex(ip);
    CEOS{k}   = ceos(ip);
    EOS{k}   = eos(ip);
    mcoef{k} = make_coeff(ip);
    Gref{k}  = Gr(ip);
end