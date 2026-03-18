function [em_comps,em_names,em_props,p_eqn] = postprocess_reactions(T,P,td,pc_id,p_out)
g0       = tl_g0(T,P,td);
em_comps = cell2mat({td(pc_id).n_em}');
em_names = cell(size(em_comps,1),1);
em_props = zeros(size(em_comps,1),1);
mu_lp    = zeros(size(em_comps,1),1);
p_eqn    = zeros(numel(pc_id),size(em_comps,1));
np       = zeros(1,numel(pc_id));
for i_p = 1:numel(pc_id)
    np(i_p) = numel(td(pc_id(i_p)).p_name);   
end
ip_end = cumsum(np);
ip_start = ip_end-np+1;
for i_p = 1:numel(pc_id)
    p_eqn(i_p,ip_start(i_p):ip_end(i_p)) = 1;   
end
cnt = 0;
for i_sol = 1:length(p_out)
    for i_p = 1:size(p_out{i_sol},1)
        if ~isempty(p_out{i_sol})
            cnt = cnt + 1;            
            em_names(p_eqn(cnt,:)==1) = td(i_sol).p_name;            
            em_props(p_eqn(cnt,:)==1) = p_out{i_sol}(i_p,:)';     
            mu_lp(p_eqn(cnt,:)==1) = tl_chemical_potential(T,P,td(i_sol),p_out{i_sol}(i_p,:),g0(i_sol));
        end
    end
end
v = null(em_comps','r');
disp_reactions(em_names,v);
chk = max(abs(v'*mu_lp));
disp(chk)

