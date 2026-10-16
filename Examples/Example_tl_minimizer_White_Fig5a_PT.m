clear,clf,addpath ../ ../Solutions ../Utilities ../EOS
solmod = 'Metapelite';
T     = linspace(450,700,20) + 273.15;
P        = linspace(0.05,1.4,21)*1e9;
Cname_oxy = {'SiO2','Al2O3','CaO','MgO','FeO','K2O','Na2O','TiO2','H2O'};
Cmol_oxy   = [74.363  9.496 0.297 3.862 7.565 3.045 0.605 0.662 50];nO         = 0.106;    
[Nsys,Cname] = Oxidemol2Elementalmol(Cname_oxy,Cmol_oxy); Nsys(strcmp(Cname,'O')) = Nsys(strcmp(Cname,'O')) + nO;Nsys = Nsys/sum(Nsys);   
sol_names    = {'Feldspar','Ilmenite','Garnet','Staurolite','Epidote','Margarite','Muscovite','Paragonite','Biotite','Chlorite','Chloritoid','Cordierite','Magnetite'};
pure_name = {'sph,tc-ds62','ru,tc-ds62','cor,tc-ds62','and,tc-ds62','ky,tc-ds62','sill,tc-ds62','q,tc-ds62','H2O,tc-ds62'};%get_pure_phases('tc-ds62',Cname);
phs_name = [sol_names,pure_name];phs_names_final = [sol_names,'Kfeldspar','Plagioclase','Melt',pure_name];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
td        = init_thermo(phs_name,Cname,'Metapelite');
for i_p = 1:length(phs_name),td(i_p).nc(:) = 3;end
p        = props_generate(td);     % generate endmember proportions
options.disp_ref  = 0;
options.disp_npc  = 0;
options.nref      = 150;
options.eps_dg    = 1e-3;
options.TPDmin    = 1;
options.fsolve    = 1;
options.solv_tol  = 2;
[T2d,P2d] = ndgrid(T,P);
asm_id = zeros(length(T)*length(P),length(phs_name));
%% Minimization refinement
dq = parallel.pool.DataQueue;afterEach(dq,@(~) printProgress());prog.total = length(T2d(:));prog.done = 0;setappdata(0,'progData',prog);
parfor iPT = 1:length(T2d(:))
    [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min,alph_fsolve{iPT},Npc_fsolve{iPT},pc_id_fsolve{iPT},p_fsolve{iPT},dg_fsolve(iPT),Nsys_chk(iPT)] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td,options);   
    send(dq,1);                 % one tick per iteration
end
fprintf('\n');                   % finish the line      
%% Postprocess
for iPT = 1:length(T2d(:))
    molm = molmass_fun(Cname);
    solv_tol = options.solv_tol;
    print_reactions = 1;
    if ~isempty(alph_fsolve{iPT}) && sum(alph_fsolve{iPT}<0)==0 && isreal(alph_fsolve{iPT})
        asm_id(iPT,1:length(pc_id_fsolve{iPT})) = pc_id_fsolve{iPT};
        Stablephases    = phs_name(pc_id_fsolve{iPT})';
        td_stable       = td(pc_id_fsolve{iPT});
        [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T2d(iPT),P2d(iPT),td_stable,alph_fsolve{iPT},Npc_fsolve{iPT},molm,p_fsolve{iPT},1:numel(pc_id_fsolve{iPT}),Stablephases,solv_tol,'PS94','S14');
        postprocess_reactions(T2d(iPT),P2d(iPT),td_stable,pc_id,p_out,print_reactions);
    else        
        [pc_id,phi,Cwt,sfu,rho,mu,p_out,phiw,g,molphase] = postprocess_fun(T2d(iPT),P2d(iPT),td,alph_all{iPT},Npc_all{iPT},molm,p_ref{iPT},pc_id_ref{iPT},phs_name,solv_tol,'PS94','S14');
        postprocess_reactions(T2d(iPT),P2d(iPT),td,pc_id,p_out,print_reactions);
        Stablephases    = phs_name(pc_id)';
        asm_id(iPT,1:length(pc_id)) = pc_id;
    end
    mineral_names = phs_name(pc_id);
    pc_id_final = zeros(size(pc_id));
    for i_p = 1:numel(pc_id)
        if numel(p_fsolve{iPT}{i_p})>1 % If solution model
            %% Renaming
            if strcmp(mineral_names(i_p),'Feldspar') % case of feldspar
                if p_fsolve{iPT}{i_p}(strcmp(td(pc_id(i_p)).p_name,'san,tc-ds62'))>0.05
                    mineral_names{i_p} = 'Kfeldspar';
                else
                    mineral_names{i_p} = 'Plagioclase';
                end
            elseif strcmp(mineral_names(i_p),'Melt(W14)') % case of melt
                mineral_names{i_p} = 'Melt';
            end
            z_all{pc_id(i_p)}(iPT,:)  = p_fsolve{i_p}*td(pc_id(i_p)).zt;
        end
        pc_id_final(i_p) = find(strcmp(phs_names_final,mineral_names(i_p)));

    end
    asm_id(iPT,1:numel(pc_id_final)) = pc_id_final;
    phs_modes(iPT,pc_id_final) = alph_ref;
    % Here some more human readable meaning of the output from the postprocess function
    Volumefractions = phi;
    Weightfractions = phiw;
    Molfractions    = molphase/sum(molphase);
    Densities       = rho; % in kg/m3
    GibbsEnergies   = g;  % in J/mol
    % % Make table
    % EquilibriumResult = [table(Volumefractions,Weightfractions,Molfractions,rho,GibbsEnergies,'RowNames',Stablephases) array2table(sfu','VariableNames',Cname)];
    % % Print in command window
    % EquilibriumResult    
    if sum(strcmp(Stablephases,'Garnet'))>0
        garnet_mode(iPT) = Volumefractions(strcmp(Stablephases,'Garnet'));
    else
        garnet_mode(iPT) = 0;
    end
    disp(iPT/length(T2d(:)))
end
figure(2)
contour(T2d-273.15,P2d/1e8,reshape(garnet_mode,length(T),length(P)),[1e-10,1e-10],'k')%,colorbar
figure(1),
tl_psection(T-273.15,P/1e9,Cname,asm_id,phs_name);
ylabel('P (GPa)'),xlabel('T(\circC)')
%% Functions
function printProgress(~)
% Retrieve and update the persistent counter
prog = getappdata(0,'progData');
prog.done = prog.done + 1;
setappdata(0,'progData',prog);       % store back

percent = (prog.done / prog.total) * 100;

% \r returns the cursor to the start of the line; fprintf overwrites it.
fprintf('\rProgress: %3.0f%% (%d/%d)', percent, prog.done, prog.total);
drawnow limitrate;   % make sure the text appears promptly
end