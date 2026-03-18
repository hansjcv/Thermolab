function [alph,Npc,pc_id,p,gmin_ref,alph_fsolve,Npc_fsolve,pc_id_fsolve,p_fsolve] = tl_minimizer(T,P,Nsys,phs_name,p,td,usr_options,rho_w,eps_di,g0,v0)
options.alph_tol    = 0;   % tolerance for alph counted as stable phase
options.nref        = 150; % max number of iterations
options.eps_dg      = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.dz_tol      = 1e-14; % tolerance to stop iterations when z window becomes below this
options.z_window    = ones(size(phs_name))*0.5; % the window over which the refined grid is generated
options.ref_fact    = 1.9; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
options.disp_ref    = 0; % show refinement graphically
options.disp_npc    = 0; % show how many pseudocompounds are generated
options.solver      = 0;
options.show_react  = 0;
options.x0          = 0.01;
options.use_pgrid   = 0;
options.gridmin     = 0;
options.x0_old      = 0;
options.algorithm   = 'dual-simplex-highs';
options.displaytype = 'off';
options.fsolve      = 0;
options.TPDmin      = 0; % use tangent plane distance and internal minimization
options.solv_tol    = 0.1;
options.fsolve_disp = 'off';
if ~exist('rho_w','var'),[rho_w,eps_di] = water_props(T,P,phs_name); end
if ~exist('eps_di','var'),[rho_w,eps_di] = water_props(T,P,phs_name); end
if exist('usr_options','var')
    option_fields = fieldnames(usr_options);
    for i = 1:length(option_fields)
        if isfield(options, option_fields{i})
            options.(option_fields{i}) = usr_options.(option_fields{i});
        else
            error('Unknown option: %s', option_fields{i});
        end
    end
end
if options.TPDmin == 1
    options.dz_tol    = 1e-200;
    if options.use_pgrid == 0
        for ip = 1:length(phs_name),p{ip} = eye(numel(td(ip).p_name));end
    end
else
    options.use_pgrid   = 1;
end
nref = options.nref;eps_dg = options.eps_dg;dz_tol=options.dz_tol;z_window = options.z_window; ref_fact = options.ref_fact; disp_ref = options.disp_ref;disp_npc = options.disp_npc;lp_algorithm=options.algorithm;disp_type = options.displaytype;
alph_tol = options.alph_tol; x0_ini = options.x0;
% Minimization refinement
if ~exist('g0','var')
    [g0,v0] = tl_g0(T,P,td,rho_w,eps_di);
end
p_it     = cell(1,length(phs_name));
td_ini   = td;
p_ini    = p;
gmin_old = 1e10;
g0_ini   = g0;
for i_sol = 1:length(phs_name)
    x0{i_sol} = ones(numel(td(i_sol).site_var),1)*x0_ini;
end
for i_ref = 1:nref    
    gmin_int = ones(size(phs_name))*1;
    [g,Npc,pc_id] = tl_gibbs_energy(T,P,phs_name,td,p,g0_ini,v0,rho_w,eps_di);    
    g = g/8.3144/T;
    LB            = zeros(1,size(g,1));
    if options.solver == 1
        Npc(Npc<1e-12)      =  0;  % Remove small number to avoid glpk instability
        [alph,gmin]  =  glpk(g,Npc,Nsys,LB,[],repmat('S',1,size(Npc,1)));
        if ~isempty(alph),exitflag(i_ref)=1;end
    elseif options.solver == 2
        [alph,gmin,exitflag(i_ref)] = opti_clp([],g,Npc,Nsys,Nsys,LB);
    else
        [alph,gmin,exitflag(i_ref),~,lambda] = linprog(g,[],[],Npc,Nsys,LB,[],optimset('Display',disp_type,'Algorithm',lp_algorithm));
    end
    if exitflag(i_ref) ~= 1,break,end
    gmin_ref(i_ref) = gmin*1e3;alph_iref{i_ref} = alph;p_iref{i_ref} = p; Nphs_iref{i_ref} = Npc; psc_id_iref{i_ref} = pc_id; g_iref{i_ref} = g;
    dg_it  = 1e8;
    dg_it1 = 1e8;
    if i_ref<nref && dg_it>eps_dg || max(z_window)>dz_tol
        p_ref = cell(1,length(phs_name));
        for i_sol = 1:length(phs_name)
            alph_ip     = alph(pc_id==i_sol);
            z{i_sol}    = p{i_sol}*td(i_sol).zt;  z{i_sol}(z{i_sol}<1+td(i_sol).z_tol & z{i_sol}>1-td(i_sol).z_tol) = 1;z{i_sol}(z{i_sol}<  td(i_sol).z_tol & z{i_sol}> -td(i_sol).z_tol) = 1e-20;
            if size(z{i_sol},1)>1
                p_it{i_sol} = [p_it{i_sol}; p{i_sol}(alph_ip>alph_tol,:)];
                if sum(alph_ip)>alph_tol
                    if options.TPDmin == 1
                        %% Adjust gibbs energy of endmembers
                        g0{i_sol} = g0_ini{i_sol} + (td(i_sol).n_em*lambda.eqlin*8.3144*T)';
                        fun1 = @(x) g_fun(x,g0(i_sol),T,P,phs_name(i_sol),td(i_sol),v0(i_sol),rho_w(i_sol),eps_di(i_sol));                        
                        g_x0 = g_fun(x0{i_sol},g0(i_sol),T,P,phs_name(i_sol),td(i_sol),v0(i_sol),rho_w(i_sol),eps_di(i_sol));
                        if ~isreal(g_x0) || isnan(g_x0) || isinf(g_x0) || options.gridmin == 1                         
                            p_grid = props_generate(td(i_sol));
                            [g_grid,~,~,~,z_grid] = tl_gibbs_energy(T,P,phs_name(i_sol),td(i_sol),p_grid,g0(i_sol),v0(i_sol),rho_w(i_sol),eps_di(i_sol));
                            [~,id] = min(g_grid);
                            z_min = z_grid(id,:);
                            x0{i_sol} = z_min(td(i_sol).site_var)';
                        end
                        try
                            x1{i_sol} = fmincon(fun1,x0{i_sol},[],[],[],[],td_ini(i_sol).z_lim(1,:),td_ini(i_sol).z_lim(2,:),[],optimoptions('fmincon','Display','off'));
                            gmin_int(i_sol) = g_fun(x1{i_sol},g0(i_sol),T,P,phs_name(i_sol),td(i_sol),v0(i_sol),rho_w(i_sol),eps_di(i_sol));
                            p_ref{i_sol} = (td(i_sol).p_from_z_cons*[1; x1{i_sol}])';                            
                            x0{i_sol} = x1{i_sol};
                        catch
                            p_ref{i_sol} = p{i_sol};
                            z_pc = z{i_sol}(alph_ip>alph_tol,:);
                            for i_pc = 1:size(z_pc,1)
                                td(i_sol).subdtype(:) = 0; % make always linear subdtype during refinement
                                td(i_sol).z_lim(1,:) = z_pc(i_pc,td(i_sol).site_var) - z_window(i_sol);
                                td(i_sol).z_lim(2,:) = z_pc(i_pc,td(i_sol).site_var) + z_window(i_sol);
                                td(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:)) = td_ini(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:));
                                td(i_sol).z_lim(1,td(i_sol).z_lim(1,:)<0) = 0;
                                p_ref{i_sol}  = cell2mat(props_generate(td(i_sol)));
                            end
                        end
                    else
                        z_pc = z{i_sol}(alph_ip>alph_tol,:);
                        for i_pc = 1:size(z_pc,1)
                            td(i_sol).subdtype(:) = 0; % make always linear subdtype during refinement
                            td(i_sol).z_lim(1,:) = z_pc(i_pc,td(i_sol).site_var) - z_window(i_sol);
                            td(i_sol).z_lim(2,:) = z_pc(i_pc,td(i_sol).site_var) + z_window(i_sol);
                            td(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:)) = td_ini(i_sol).z_lim(2,td(i_sol).z_lim(2,:)>td_ini(i_sol).z_lim(2,:));
                            td(i_sol).z_lim(1,td(i_sol).z_lim(1,:)<0) = 0;
                            p_ref{i_sol}  = [p_ref{i_sol}; cell2mat(props_generate(td(i_sol)))];
                            % p_ref{i_sol}  = cell2mat(props_generate(td(i_sol)));
                        end
                    end
                else
                    p_ref{i_sol} = p{i_sol}; % this is in some cases important as it keeps also the unstable phases
                end
            else
                x1{i_sol} = 1;
                if options.TPDmin ~= 1                    
                    p_ref{i_sol} = p{i_sol};
                end
            end
        end
    end
    z_window = z_window/ref_fact;
    if i_ref>3
        dg_it = max(abs((diff(gmin_ref(i_ref-3:i_ref)))));
        dg_it1 = diff(gmin_ref(i_ref-1:i_ref));
    end
    if options.show_react == 1
        [~,~,p_out,pc_id_out] = cluster_p(alph,Npc,p,pc_id,1,phs_name);
        postprocess_reactions(T,P,td,pc_id_out,p_out,[]);
    end
    if i_ref>3 && dg_it<eps_dg || max(z_window)<dz_tol,break,end
    if i_ref < nref && gmin_old >= gmin
        for i_sol = 1:length(phs_name)
            if options.TPDmin == 1
                p{i_sol} = [p_ini{i_sol};p_it{i_sol}; p_ref{i_sol}];                
            else
                p{i_sol} = [p_it{i_sol}; p_ref{i_sol}];
            end
            if disp_npc == 1
                disp([phs_name{i_sol} ' ' num2str(size(p{i_sol},1))]);
            end
        end
    end
    if disp_ref == 1
        if options.TPDmin == 1
            plot(1:i_ref,gmin_ref(1:i_ref)),title([dg_it1,i_ref,length(pc_id)]),drawnow
        else
            plot(1:i_ref,gmin_ref(1:i_ref)),title([dg_it,i_ref,length(pc_id),max(z_window)]),drawnow%hold on,text(i_ref,gmin_ref(i_ref),phs_name(pc_id)),            
        end
    end
    gmin_old = gmin;
end
% Postprocess
alph     = alph_iref{find(exitflag==1,1,'last')};
p        = p_iref{find(exitflag==1,1,'last')};
Npc      = Nphs_iref{find(exitflag==1,1,'last')};
pc_id    = psc_id_iref{find(exitflag==1,1,'last')};
gmin_ref = gmin_ref(1:find(exitflag==1,1,'last'));
for ip = 1:length(p)
    unstb_id = alph(pc_id==ip)<=alph_tol;
    p{ip}(unstb_id,:) = []; % throw out zeros
end
solv_tol = options.solv_tol;
Npc(:,alph<=alph_tol)   = []; % throw out zeros
pc_id(alph<=alph_tol)   = []; % throw out zeros
alph(alph<=alph_tol)    = []; % throw out zeros
if options.fsolve == 1    
    [alph_out,Npc_out,p_out,pc_id_out] = cluster_p(alph,Npc,p,pc_id,solv_tol,phs_name);
    p_fsolve          = cell(size(pc_id_out));
    pc_id_fsolve      = pc_id_out;
    [em_comps,em_names,em_props,p_eqn] = postprocess_reactions(T,P,td_ini,pc_id_out,p_out);
    v = null(em_comps','r');
    options_fsolve = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display',options.fsolve_disp);%'MaxFunctionEvaluations',10000,'MaxIterations',10000,'StepTolerance',1e-16,'TolFun',1e-16
    x0_fsolve      = [em_props(1:sum(sum(p_eqn(sum(p_eqn,2)>1,:))));alph_out];
    f              = @(x) tl_dmu(x,T,P,phs_name(pc_id_out),Nsys,td_ini(pc_id_out),p_eqn,g0_ini(pc_id_out),v,em_comps);
    x              = fsolve(f,x0_fsolve,options_fsolve); % find proportions
    alph_fsolve       = x(sum(sum(p_eqn(sum(p_eqn,2)>1,:)))+1:sum(sum(p_eqn(sum(p_eqn,2)>1,:))) + numel(pc_id_out));
    p_ref          = x(1:sum(sum(p_eqn(sum(p_eqn,2)>1,:))));
    %% Postprocessing 2
    em_props   = [p_ref;ones(sum(sum(p_eqn,2)==1),1)];
    em_chempot = zeros(size(em_props));
    Npc_fsolve = zeros(numel(Nsys),numel(pc_id_out));
    for i_p = 1:length(pc_id_out)
        p_fsolve{i_p}    = em_props(p_eqn(i_p,:)==1)';
        Npc_fsolve(:,i_p)           = p_fsolve{i_p}*em_comps(p_eqn(i_p,:)==1,:);
        em_chempot(p_eqn(i_p,:)==1) = tl_chemical_potential(T,P,td_ini(pc_id_out(i_p)),p_fsolve{i_p},g0_ini(pc_id_out(i_p)))';
    end    
    chk_dg_ref   = max(abs(v'*em_chempot))
    chk_Nsys_ref = max(abs(Npc_fsolve*alph_fsolve-Nsys'))    
else
    alph_fsolve  = [];
    p_fsolve     = [];
    Npc_fsolve   = [];
    pc_id_fsolve = [];
end
if options.TPDmin == 1
    plot(1:length(gmin_ref),gmin_ref),title([dg_it1,i_ref,length(pc_id)]),drawnow
else
    plot(1:length(gmin_ref),gmin_ref),title([dg_it,i_ref,length(pc_id),max(z_window)]),drawnow
end
end
%% Functions
function g = g_fun(x,g0,T,P,phs_name,td,v0,rho_w,eps_di)
p{1}    = (td.p_from_z_cons*[1; x])';
g       = tl_gibbs_energy(T,P,phs_name,td,p,g0,v0,rho_w,eps_di)/8.3144/T;
end
function Feqn = tl_dmu(x,T,P,mineral_names,Nsys,td_eqn,p_eqn,g0_eqn,v,em_comps)
em_chempot = zeros(size(p_eqn,2),1);
em_props = [x(1:sum(sum(p_eqn(sum(p_eqn,2)>1,:))));ones(sum(sum(p_eqn,2)==1),1)];
alph_out = x(sum(sum(p_eqn(sum(p_eqn,2)>1,:)))+1:sum(sum(p_eqn(sum(p_eqn,2)>1,:)))+length(mineral_names));
for i_p = 1:length(mineral_names)    
    em_chempot(p_eqn(i_p,:)==1) = tl_chemical_potential(T,P,td_eqn(i_p),em_props(p_eqn(i_p,:)==1)',g0_eqn(i_p))';
end
Feqn = [p_eqn(sum(p_eqn,2)>1,1:sum(sum(p_eqn(sum(p_eqn,2)>1,:))))*em_props(1:sum(sum(p_eqn(sum(p_eqn,2)>1,:)))) - ones(size(p_eqn(sum(p_eqn,2)>1,:),1),1);
        v'*em_chempot;
        (p_eqn*(em_props.*em_comps))'*alph_out-Nsys'];
end