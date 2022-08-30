function [alph_level2, T_level2,P_level2,refine_id,Npc_level2,pc_level2,p_level2] = refine_grid(T,P,alph,asm_id,p,pc_id,Npc)

% alph = reshape(alph(:),length(T)*length(P),length(psc_id));

% Grid stuff
nx = length(T);
ny = length(P);

% Node to global (like FEM, square linear elements)
n2g = zeros(4,(nx-1)*(ny-1));
for j = 1:ny-1
    for i = 1:nx-1
        iel        = i   + (j-1)*(nx-1);
        n2g(1,iel) = iel + (j-1);
        n2g(2,iel) = n2g(1,iel)+1;
        n2g(3,iel) = n2g(2,iel)+nx;
        n2g(4,iel) = n2g(3,iel)-1;
    end
end
% Nodes to level1 elements
n2el2 = zeros(4,(nx-1)*(ny-1));
for j = 1:ny-1
    for i = 1:nx-1
        iel   = i   + (j-1)*(nx-1);
        n2el2(1,iel) = 2*i - 1 + (nx-1)*(j-1)*4;
        n2el2(2,iel) = n2el2(1,iel) + 1;
        n2el2(3,iel) = n2el2(1,iel) + 2*(nx-1);
        n2el2(4,iel) = n2el2(3,iel) + 1;
    end
end
% Subdivide grid
nx = nx*2-1;
ny = ny*2-1;
dT = (max(T)-min(T))/(nx-1);
dP = (max(P)-min(P))/(ny-1);
T2 = min(T):dT:max(T);
P2 = min(P):dP:max(P);
[T2d2,P2d2 ]= ndgrid(T2,P2);
% Nodes to global level2
n2g2 = zeros(4,(nx-1)*(ny-1));
for j = 1:ny-1
    for i = 1:nx-1
        iel        = i   + (j-1)*(nx-1);
        n2g2(1,iel) = iel + (j-1);
        n2g2(2,iel) = n2g2(1,iel)+1;
        n2g2(3,iel) = n2g2(2,iel)+nx;
        n2g2(4,iel) = n2g2(3,iel)-1;
    end
end
% Check assemblage at each level1 element
diff_ass = zeros(size(n2g'));
for i = 1:size(n2g,1)
    diff_ass(:,i) = asm_id(n2g(1,:))-asm_id(n2g(i,:));
end
chk1       =  (sum(abs(diff_ass),2)==0);
elem_done   = find(chk1); % cells that are completely finished
node_done   = n2g(:,elem_done); % nodes that are completely finished
elem2_done  = n2el2(:,chk1); % new smaller cells that are completely finished
% Save the finished assemblages in new array
phi_all2       = cell(length(T2)*length(P2),1);
pc_all2        = cell(length(T2)*length(P2),1);
Npc_all2        = cell(length(T2)*length(P2),1);
p_all2        = cell(length(T2)*length(P2),1);
phi_elem       = cell(length(elem_done),1);
pc_elem       = cell(length(elem_done),1);
Npc_elem       = cell(length(elem_done),1);
p_elem       = cell(length(elem_done),1);
for i = 1:length(elem_done)
    phi_elem(i) = alph(node_done(1,i));
    pc_elem(i)  = pc_id(node_done(1,i));
    p_elem(i)   = p(node_done(1,i));
    Npc_elem(i) = Npc(node_done(1,i));
    nodes_i     = n2g2(:,n2el2(:,elem_done(i))); 
%     plot(T2d(n2g(:,elem_done)),P2d(n2g(:,elem_done)),'w*'),hold on    
    for j = 1:length(nodes_i(:))
        phi_all2(nodes_i(j)) = phi_elem(i);
        pc_all2(nodes_i(j))  = pc_elem(i);
        Npc_all2(nodes_i(j)) = Npc_elem(i);
        p_all2(nodes_i(j))   = p_elem(i);
%         plot(T2d2(nodes_i(j)),P2d2(nodes_i(j)),'wo')
    end
end

% Save the finished assemblages in new array
% phi_all2       = zeros(length(T2)*length(P2),length(psc_id));
% phi_elem       = zeros(length(elem_done),length(psc_id));
% for i = 1:length(elem_done)
%     phi_elem(i,:) = alph(node_done(1,i),:);
%     nodes_i       = n2g2(:,n2el2(:,elem_done(i))); 
% %     plot(T2d(n2g(:,elem_done)),P2d(n2g(:,elem_done)),'w*'),hold on    
%     for j = 1:length(nodes_i(:))
%         phi_all2(nodes_i(j),:) = phi_elem(i,:);
% %         plot(T2d2(nodes_i(j)),P2d2(nodes_i(j)),'wo')
%     end
% end


T2d2_ind = ones(size(T2d2(:)));
T2d2_ind(n2g2(:,elem2_done(:))) = 0;
T2d2_ind = logical(T2d2_ind);

T_level2    = T2;
P_level2    = P2;
alph_level2 = phi_all2;
refine_id   = T2d2_ind;
Npc_level2  = Npc_all2;
pc_level2  = pc_all2;
p_level2  = p_all2;
% hold off
