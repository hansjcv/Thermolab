function [alph,Npc,pc_id,p,gmin_ref] = tl_minimizer(T,P,Nsys,phs_name,p,td,options,rho_w,eps_di,g0,v0)
alph_tol = 0;   % tolerance for alph counted as stable phase
if ~exist('rho_w','var'),[rho_w,eps_di] = water_props(T,P,phs_name); end
if ~exist('eps_di','var'),[rho_w,eps_di] = water_props(T,P,phs_name); end
if ~exist('options','var')
    options.nref     = 150; % max number of iterations
    options.eps_dg   = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
    options.dz_tol   = 1e-14; % tolerance to stop iterations when z window becomes below this
    options.z_window = ones(size(phs_name))*0.085; % the window over which the refined grid is generated
    options.dz_fact  = ones(size(phs_name))*1.5; % the factor to determine new dz spacing, the larger, the more pseudocompounds
    options.ref_fact = 1.25; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
    options.disp_ref = 0; % show refinement graphically
    options.solver   = 0;
end
nref = options.nref;eps_dg = options.eps_dg;dz_tol=options.dz_tol;z_window = options.z_window; dz_fact = options.dz_fact; ref_fact = options.ref_fact; disp_ref = options.disp_ref;
% Minimization refinement
if ~exist('g0','var')
    [g0,v0] = tl_g0(T,P,td,rho_w,eps_di);
end
p_it    = cell(1,length(phs_name));
td_ini  = td;
gmin_old = 1e10;
for i_ref = 1:nref       
    [g,Npc,pc_id] = tl_gibbs_energy(T,P,phs_name,td,p,g0,v0,rho_w,eps_di);    
    LB            = zeros(1,size(g,1));    
    Npc(Npc<1e-14&Npc>-1e-14)      =  0;
    if options.solver == 1
        Npc(Npc<1e-12)      =  0;  % Remove small number to avoid glpk instability
        [alph,gmin]  =  glpk(g,Npc,Nsys,LB,[],repmat('S',1,size(Npc,1)));     
        if ~isempty(alph),exitflag(i_ref)=1;end
    elseif options.solver == 2
        [alph,gmin,exitflag(i_ref)] = opti_clp([],g,Npc,Nsys,Nsys,LB);    
    else
        % g = g-g(pc_id==find(strcmp(phs_name,'q,tc-ds55')),:).*Npc(1,:)';
        % [alph,gmin,exitflag(i_ref)] = linprog(g,[],[],Npc(2:end,:),Nsys(2:end),LB,[],optimset('Display','off'));            
        [alph,gmin,exitflag(i_ref)] = linprog(g,[],[],Npc,Nsys,LB,[],optimset('Display','off'));    
    end    
    if exitflag(i_ref) ~= 1,break,end
    gmin_ref(i_ref) = gmin;alph_iref{i_ref} = alph;p_iref{i_ref} = p; Nphs_iref{i_ref} = Npc; psc_id_iref{i_ref} = pc_id; g_iref{i_ref} = g;    
    dg_it  = 1e8;       
    if i_ref<nref && dg_it>eps_dg || max(z_window)>dz_tol
        p_ref = cell(1,length(phs_name));
        % this is to do clustering
%         for i_dz = 1:length(phs_name)
%             if ~isempty(td(i_dz).dz)
%                 solv_tol(i_dz) = max(td(i_dz).dz)*1.5;
%             end
%         end
%         solv_tol = max(solv_tol);
%         [alph_out,Npc_out,p_out,pc_id_out] = cluster_p(alph,Npc,p,pc_id,solv_tol,phs_name);
        %
        for i_sol = 1:length(phs_name)
            alph_ip     = alph(pc_id==i_sol);            
            z{i_sol}    = p{i_sol}*td(i_sol).zt;  z{i_sol}(z{i_sol}<1+td(i_sol).z_tol & z{i_sol}>1-td(i_sol).z_tol) = 1;z{i_sol}(z{i_sol}<  td(i_sol).z_tol & z{i_sol}> -td(i_sol).z_tol) = 1e-20;
            if size(z{i_sol},1)>1
                p_it{i_sol} = [p_it{i_sol}; p{i_sol}(alph_ip>0,:)];
                if sum(alph_ip)>alph_tol
                    z_pc = z{i_sol}(alph_ip>alph_tol,:);    
%                     z_pc = p_out{i_sol}*td(i_sol).zt;
                    for i_pc = 1:size(z_pc,1)
                        td(i_sol).subdtype(:) = 0;
                        td(i_sol).z_lim(1,:) = z_pc(i_pc,td(i_sol).site_var) - z_window(i_sol);
                        td(i_sol).z_lim(2,:) = z_pc(i_pc,td(i_sol).site_var) + z_window(i_sol);
                        td(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:)) = td_ini(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:));
                        td(i_sol).z_lim(1,td(i_sol).z_lim(1,:)<0) = 0;                        
                        td(i_sol).dz  = abs(diff(td(i_sol).z_lim))./dz_fact(i_sol);
                        p_ref{i_sol}  = [p_ref{i_sol}; cell2mat(props_generate(td(i_sol)))];
                    end
                else
                    p_ref{i_sol} = p{i_sol}; % this is in some cases important as it keeps also the unstable phases
                end
            else
                p_ref{i_sol} = p{i_sol};
            end
        end
        z_window = z_window/ref_fact;        
        if i_ref>3
            dg_it = max(abs((diff(gmin_ref(i_ref-3:i_ref)))));
        end        
        if i_ref>3 && dg_it<eps_dg || max(z_window)<dz_tol,break,end        
        if i_ref < nref && gmin_old >= gmin
            for i_sol = 1:length(phs_name)
                p{i_sol} = [p_it{i_sol}; p_ref{i_sol}];
            end
        end        
    end    
    if disp_ref == 1
        plot(1:i_ref,gmin_ref(1:i_ref)),title([dg_it,i_ref,length(pc_id),max(z_window)]),drawnow
    end           
    gmin_old = gmin;
end
alph = alph_iref{find(exitflag==1,1,'last')};
p = p_iref{find(exitflag==1,1,'last')};
Npc = Nphs_iref{find(exitflag==1,1,'last')};
pc_id = psc_id_iref{find(exitflag==1,1,'last')};
gmin_ref = gmin_ref(1:find(exitflag==1,1,'last'));
for ip = 1:length(p)
    unstb_id = alph(pc_id==ip)==0;
    p{ip}(unstb_id,:) = []; % throw out zeros    
end
Npc(:,alph==0)   = []; % throw out zeros
pc_id(alph==0)   = []; % throw out zeros
alph(alph==0)    = []; % throw out zeros
plot(1:length(gmin_ref),gmin_ref),title([dg_it,i_ref,length(pc_id),max(z_window)]),drawnow