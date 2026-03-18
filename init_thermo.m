function td = init_thermo(mineral,Cname,solfile,excel_yes,pure_names,dz_fluid,max_frac_fl)
load 'tl_dataset' nphs elements phs_names make_coeff td_data ceos eos dgex Gr rho_w_mod dielec_mod dHf
% if ~exist(pure_names)
% end
names_table = readtable('mineral_list.xlsx');
if ~exist('solfile','var')
    sol_mod_list = [];
else
    if ~exist('excel_yes','var')
        excel_yes = 1;
    end
    if excel_yes == 1
        if exist('sheetnames','builtin')==2
            sol_mod_list = sheetnames([solfile '.xlsx']);
        else
            [~,sol_mod_list] = xlsfinfo([solfile '.xlsx']);
        end
    else
        load(solfile,'data_all','txt_all','sol_mod_list');
    end
end
if exist('Cname','var')
    c_ind = find_c_ind(Cname,elements);
    oxide_flag = 0;
    if sum(c_ind>0)==0
        Cname_oxide = Cname;
        [~,Cname] = Oxidemol2Elementalmol(Cname_oxide,zeros(size(Cname_oxide)));
        c_ind = find_c_ind(Cname,elements);
        oxide_flag = 1;
    end
else
    c_ind = 1:length(elements);
    oxide_flag = 0;
end
for i_sol = 1:length(mineral)
    if ~strcmp(mineral(i_sol),'Aq-fluid')
        Nphs = []; chg =[];p_id = []; em_data = [];dGex = []; CEOS = []; EOS = []; mcoef = []; Gref = [];w=[];alp=[];
        if sum(strcmp(sol_mod_list,mineral{i_sol}))>0
            % Read data from excel tables
            if excel_yes == 1
                [data,txt]  = xlsread([solfile '.xlsx'],mineral{i_sol},'','basic');
            else
                data = data_all{strcmp(sol_mod_list,mineral{i_sol})};
                txt  = txt_all{strcmp(sol_mod_list,mineral{i_sol})};
            end
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
            nc              = (z_max-z_min)./(dz);%ones(size(dz))*3;
            nc = fix((z_max-z_min)./dz+1);
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
            p_name = {mineral{i_sol}};
            alp = [
                1.000
                0.000
                0.000
                ];
            w0      = 0;
            wT      = 0;
            wP      = 0;
            st      = 1;
            site_id = 1;
            mtpl    = 1;
            mod_id = 0;
            z_lim      = [];
            dz         = 1;
            subdtype  = 1;
            nc = 1;
        end
        for k = 1:length(p_name)
            ip  = find(strcmp(phs_names,p_name(k)));                 % Find index of the endmember in the list
            if isempty(ip)                 
                phs_abbr = names_table.Abbreviation{strcmp(names_table.Mineral,p_name(k))};
                ip  = find(strcmp(phs_names,[phs_abbr,',tc-ds633'])); 
            end
            Nphs(k,:) = make_coeff(ip)*nphs(ip,:);
            chg(k)    = make_coeff(ip)*nphs(ip,strcmp(elements,'e'));
            p_id{k} = ip;
            em_data{k}   = cell2mat(td_data(ip)');
            dGex{k}   = dgex(ip);
            CEOS{k}   = ceos(ip);
            EOS{k}   = eos(ip);
            mcoef{k} = make_coeff(ip);
            Gref{k}  = Gr(ip);
            dHform{k} = dHf(ip);
            rho_w_m{k} = rho_w_mod(ip);
            dielec_w_m{k} = dielec_mod(ip);
        end
        if strcmp(mineral{i_sol},'Melt(H18)')
            z_name = {'pq','psl','pwo','pjd','phm','pek','pti','pkj','pct','pol','sumT','mgM','feM','CaM','AlM','sumM','xh','xv'};
            z_fac = {'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)'  '4./(3*yct+4)' '4./(3*yct+4)' '4./(3*yct+4)' '4./(3*yct+4)' '4./(3*yct+4)' '4./(3*yct+4)' '4./(3*yct+4)'};
            mtpl = {'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)' 'ones(size(p,1),1)'  'ones(size(p,1),1)' '-1*ones(size(p,1),1)' '((3*yct)/4 + 1)', '((3*yct)/4 + 1)', '((3*yct)/4 + 1) - yct./CaM', '((3*yct)/4 + 1) - yct./AlM', '-((3*yct)/4 + 1) + 2*yct./sumM', '(3*yct + 4)/2', '((3*yct)/2 + 2)-(3/2*yct)./xv'};
            for iz = 1:length(occupancy)
                z_name_id(iz) = find(strcmp(z_name,occupancy(iz)));
            end
            z_name = z_name(z_name_id);
            z_fac = z_fac(z_name_id);
            mtpl = mtpl(z_name_id);
        else
            z_name = [];
            z_fac  = [];
        end
        if sum(c_ind==-1)>0
            inc_phs_id  = sum(Nphs(:,c_ind==-1),2)==0;
            p_id   = p_id(inc_phs_id);
            p_name = p_name(inc_phs_id);
            alp    = alp(:,inc_phs_id);
            w0     = w0(inc_phs_id,:,:);
            w0     = w0(:,inc_phs_id,:);
            wT     = wT(inc_phs_id,:,:);
            wT     = wT(:,inc_phs_id,:);
            wP     = wP(inc_phs_id,:,:);
            wP     = wP(:,inc_phs_id,:);
            st     = st(inc_phs_id,:);
            site_id(sum(abs(st))==0) = [];
            occupancy(sum(abs(st))==0) = [];
            sitenames(sum(abs(st))==0) = [];
            z_lim(:,sum(abs(st))==0) = [];
            mtpl(sum(abs(st))==0) = [];
            if strcmp(mineral{i_sol},'Melt(H18)')
                z_fac(sum(abs(st))==0)  = [];
                z_name(sum(abs(st))==0) = [];
            end
            st(:,sum(abs(st))==0) = [];
            Nphs(~inc_phs_id,:) = [];
            chg(~inc_phs_id) = [];
            em_data(~inc_phs_id) = [];
            dGex(~inc_phs_id)   = [];
            CEOS(~inc_phs_id)   = [];
            EOS(~inc_phs_id)   = [];
            mcoef(~inc_phs_id) = [];
            Gref(~inc_phs_id)  = [];
            dHform(~inc_phs_id)  = [];
            rho_w_m(~inc_phs_id) = [];
            dielec_w_m(~inc_phs_id) = [];
        end
        n_em = zeros(size(Nphs,1),sum(c_ind~=-1));
        n_em(:,c_ind(c_ind~=-1))    = Nphs(:,c_ind~=-1);
        if oxide_flag == 1
            n_em_ox= zeros(size(Nphs,1),numel(Cname_oxide));
            for i_em = 1:size(n_em,1)
                n_em_ox(i_em,:) = Elementalmol2Oxidemol(Cname,n_em(i_em,:));
            end
            n_em = n_em_ox;
        end
        w(:,:,1) = w0;
        w(:,:,2) = wT;
        w(:,:,3) = wP;
        [zt,zmax]                               = st2zt(st,site_id);    % Site fraction table
        [p_from_c_cons,icomp_indep,isite_indep] = comp2prop(n_em',zt'); % C variables, p to C conversion matrix
        [p_from_z_cons,site_var,cons] = sites2prop(zt);
        if mod_id == 4
            [p_from_c_cons,icomp_indep,isite_indep] = comp2prop(n_em',st'); % C variables, p to C conversion matrix
            [p_from_z_cons,site_var,cons] = sites2prop(st);
        end
        
        if strcmp(mineral{i_sol},'Melt(H18)')
            zt = st;
            z_fac_string = 'z_fac = [';
            mtpl_string = 'mtpl = [';
            z_name_string = [];
            for i_zfac = 1:length(z_fac)
                z_fac_string = [ z_fac_string ' ' z_fac{i_zfac} ];
                mtpl_string = [ mtpl_string ' ' mtpl{i_zfac} ];
                z_name_string = [z_name_string z_name{i_zfac} ' = z(:,' num2str(i_zfac) ');' ];
            end
            z_fac_string = [ z_fac_string '];'];
            mtpl_string = [ mtpl_string '];'];
            mtpl = mtpl_string;
            z_fac = z_fac_string;
            z_name = z_name_string;
        else
            z_fac = [];
            z_name = [];
        end
        td(i_sol).p_id    = p_id;
        td(i_sol).st      = st;
        td(i_sol).site_id = site_id;
        td(i_sol).mtpl    = mtpl;
        td(i_sol).mod_id  = mod_id;
        td(i_sol).alp     = alp;
        td(i_sol).w       = w;
        td(i_sol).n_em    = n_em;
        td(i_sol).chg     = chg;
        td(i_sol).em_data = em_data;
        td(i_sol).dGex    = dGex;
        td(i_sol).EOS     = EOS;
        td(i_sol).CEOS    = CEOS;
        td(i_sol).mcoef   = mcoef;
        td(i_sol).Gref    = Gref;
        td(i_sol).dHf     = dHform;
        td(i_sol).zt      = zt;
        td(i_sol).zmax    = zmax;
        td(i_sol).p_from_c_cons = p_from_c_cons;
        td(i_sol).icomp_indep = icomp_indep;
        td(i_sol).isite_indep = isite_indep;
        td(i_sol).p_from_z_cons = p_from_z_cons;
        td(i_sol).site_var = site_var;
        td(i_sol).cons = cons;
        td(i_sol).z_lim  = z_lim(:,site_var);
        td(i_sol).z_lims = z_lim;
        td(i_sol).dz     = dz(site_var);
        td(i_sol).nc     = nc(site_var);
        td(i_sol).subdtype     = subdtype;
        td(i_sol).phase_name    = mineral{i_sol};
        td(i_sol).z_tol         = 1e-10;
        td(i_sol).p_name        = p_name;
        td(i_sol).z_fac         = z_fac;
        td(i_sol).z_name        = z_name;
        td(i_sol).occupancy     = occupancy;
        td(i_sol).sitenames        = sitenames;
    else
        % c_ind = find_c_ind(Cname,elements);
        p_id = zeros(1,numel(pure_names));
        for ip = 1:length(pure_names)
            p_id(ip)          = find(strcmp(phs_names,pure_names{ip}));
            td(i_sol).p_id{ip}    = p_id(ip);
            td(i_sol).em_data{ip} = td_data{p_id(ip)};
            td(i_sol).dGex{ip}    = dgex(p_id(ip));
            td(i_sol).EOS{ip}     = eos(p_id(ip));
            td(i_sol).CEOS{ip}    = ceos(p_id(ip));
            td(i_sol).mcoef{ip}   = make_coeff(p_id(ip));
            td(i_sol).Gref{ip}    = Gr(p_id(ip));
            td(i_sol).dHf{ip}     = dHf(p_id(ip));
        end
        Nphs                     = nphs(p_id,:);
        n_em                     = zeros(size(Nphs,1),sum(c_ind~=-1));
        n_em(:,c_ind(c_ind~=-1)) = Nphs(:,c_ind~=-1);
        st                       = eye(numel(p_id));
        site_id                  = ones(1,numel(p_id));
        z_lim                    = [zeros(1,numel(p_id));ones(1,numel(p_id))*max_frac_fl];
        dz                       = ones(1,numel(p_id))*dz_fluid;
        nc                       = ones(1,numel(p_id))*3;
        subd_type                = ones(1,numel(p_id));
        [zt,zmax]                               = st2zt(st,site_id);    % Site fraction table
        [p_from_c_cons,icomp_indep,isite_indep] = comp2prop(n_em',zt'); % C variables, p to C conversion matrix
        [p_from_z_cons,site_var,cons]           = sites2prop(zt);
        td(i_sol).st                 = st;
        td(i_sol).site_id            = site_id;
        td(i_sol).mtpl               = ones(1,numel(p_id));
        td(i_sol).mod_id             = 5;
        td(i_sol).alp                = [ ones(1,numel(p_id));
            zeros(1,numel(p_id));
            zeros(1,numel(p_id))];
        td(i_sol).w                  = zeros(numel(p_id),numel(p_id),3);
        td(i_sol).n_em               = n_em;
        td(i_sol).chg                = n_em(:,end);
        td(i_sol).zt                 = zt;
        td(i_sol).zmax               = zmax;
        td(i_sol).p_from_c_cons      = p_from_c_cons;
        td(i_sol).icomp_indep        = icomp_indep;
        td(i_sol).isite_indep        = isite_indep;
        td(i_sol).p_from_z_cons      = p_from_z_cons;
        td(i_sol).site_var           = site_var;
        td(i_sol).cons               = cons;
        td(i_sol).z_lim              = z_lim(:,site_var);
        td(i_sol).z_lims             = z_lim;
        td(i_sol).dz                 = dz;
        td(i_sol).nc                 = nc;
        td(i_sol).subdtype           = subd_type;
        td(i_sol).phase_name         = mineral{i_sol};
        td(i_sol).z_tol              = 1.0000e-14;
        td(i_sol).p_name             = pure_names;
        td(i_sol).z_fac              = [];
        td(i_sol).z_name             = [];
        td(i_sol).occupancy          = pure_names;
        td(i_sol).sitenames          = repmat({'A'},1,numel(p_id));
    end
    if sum(abs(td(i_sol).n_em))==0,error([td(i_sol).phase_name ' not known, check for typo']),end
end