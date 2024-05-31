function phases = get_pure_phases(db_name,Cname,is_liq,is_gas,is_aq)
load('tl_dataset','elements','db_names','phs_names','nphs')
if nargin < 3
    is_aq   = 0;
    is_gas  = 0;
    is_liq  = 0;
elseif nargin < 4
    is_aq   = 0;
    is_gas  = 0;
elseif nargin <5
    is_aq = 0;
end
db_ind = strcmp(db_names,db_name);
phases = phs_names(db_ind);
nphs   = nphs(db_ind,:);
exc_id = [];
switch db_name
    case 'tc-ds55'
        phases = phases(1:190);
        nphs   = nphs(1:190,:);        
        gas_id = 151:156;
        liq_id = 157:168;
        aq_id  = 169:190;
    case 'tc-ds62'
        phases = phases(1:256);
        nphs   = nphs(1:256,:);
        liq_id = 213:230;
        gas_id = 205:212;
        aq_id  = 231:256;                
    case 'tc-ds633'
        phases = phases(1:289);
        nphs   = nphs(1:289,:);        
        gas_id = 234:241;
        liq_id = 242:263;
        aq_id  = 264:289;
    case 'DEW'
        aq_id  = 1:length(phases);
        liq_id = [];
        gas_id = [];
    case 'supcrt'
        aq_id  = 1:length(phases);
        liq_id = [];
        gas_id = [];
    case 'Miron'
        aq_id  = 1:length(phases);
        liq_id = [];
        gas_id = [];
end
if is_aq == 0, exc_id = [exc_id, aq_id]; end
if is_gas == 0,exc_id = [exc_id, gas_id]; end
if is_liq == 0,exc_id = [exc_id, liq_id]; end
phases(exc_id) = [];
nphs(exc_id,:) = [];
for i_c = 1:length(Cname)
    c_ind(i_c) = find(strcmp(elements,Cname(i_c)));
end
c_excl        = ones(1,length(elements));
c_excl(c_ind) = 0;
phases = phases(sum(abs(nphs(:,c_excl==1)),2)==0&sum(abs(nphs(:,c_excl==0)),2)>0);