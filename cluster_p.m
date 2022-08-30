function [alph,Npc,p_out,pc_id] = cluster_p(alph_in,Npc_in,p_in,pc_id_in,solv_tol,phs_id)
for ip = 1:length(phs_id)
    alph  = alph_in(pc_id_in==ip);
    Npc   = Npc_in(:,pc_id_in==ip);
    p     = p_in{ip};
    pc_id = pc_id_in(pc_id_in==ip);
    out  = 0;
    while out == 0
        i_stab = find(alph>0);
%         X_stab = Npc(:,alph>0);
        X_stab = p(alph>0,:)';
        nst    = size(X_stab,2);
        if nst > 1
            d_curr = 10000;
            for i = 1:nst
                for j = 1:nst
                    if i~=j
                        d_curr  = X_stab(:,i) - X_stab(:,j);
                        d_curr  = sqrt(sum(d_curr.*d_curr));
                    end
                    if d_curr < solv_tol,break,end
                end
                if d_curr < solv_tol,break,end
            end
            if d_curr < solv_tol
                alph_nei = alph(i_stab([i,j]));
                sum_alph = sum(alph_nei);
                alph_nei = alph_nei/sum_alph;
                Npc(:,i_stab(i)) = Npc(:,i_stab([i,j]))*alph_nei;
                p(i_stab(i),:)   = alph_nei'*p(i_stab([i,j]),:);
                alph(i_stab(i))  = sum_alph;
                alph(i_stab(j))  = 0;
            else
                out = 1;
            end
        else
            out = 1;
        end
    end
    Npc_out{ip}   = Npc(:,alph>0);
    p_out{ip}     = p(alph>0,:);%(alph>0,:);
    pc_id_out{ip} = pc_id(alph>0);
    alph_out{ip}  = alph(alph>0);
end
Npc = cell2mat(Npc_out);
pc_id = cell2mat(pc_id_out);
alph = cell2mat(alph_out');