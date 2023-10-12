function [g0,v0,s0] = tl_g0(T,P,td,rho_w,eps_di,Gapp)
if ~exist('rho_w','var'),[rho_w,eps_di] = water_props(T,P,'H2O,tc-ds55');end     
if ~exist('eps_di','var'),[rho_w,eps_di] = water_props(T,P,'H2O,tc-ds55'); end
if ~exist('Gapp','var'),Gapp = 0;end
for i_sol = 1:numel(td)
    em_data = td(i_sol).em_data;CEOS = td(i_sol).CEOS;EOS = td(i_sol).EOS;dGex = td(i_sol).dGex;Gr = td(i_sol).Gref;mcoef = td(i_sol).mcoef;dHf = td(i_sol).dHf;
    for k = 1:length(em_data)
        Gjp = zeros(numel(T),size(em_data{k},1)); % Initialize G of make endmembers
        Vjp = zeros(numel(T),size(em_data{k},1)); % Initialize G of make endmembers
        Sjp = zeros(numel(T),size(em_data{k},1)); % Initialize G of make endmembers
        Hfjp = zeros(numel(T),size(em_data{k},1)); % Initialize G of make endmembers
        for jp = 1:size(em_data{k},1)      % For endmembers made of multiple phases
            [SdT,S]   = intSdT(T,P,em_data{k}(jp,:),CEOS{k}(jp));              % SdT integral
            [VdP,V]   = intVdP(T,P,em_data{k}(jp,:), EOS{k}(jp));              % VdP integral
            Gexc      = g0_exc(T,P,em_data{k}(jp,:),dGex{k}(jp),rho_w(:,i_sol),eps_di(:,i_sol)); % dG of any fitting excess energy
            Gjp(:,jp) = Gr{k}(jp) - SdT + VdP + Gexc;                     % Equation 4
            Vjp(:,jp) = V;
            Sjp(:,jp) = S; 
            Hfjp(:,jp) = dHf{k}(jp);
        end
        g0{i_sol}(:,k)     = mcoef{k}*Gjp';                          % The Gibbs energy of the phase
        v0{i_sol}(:,k)     = mcoef{k}*Vjp';
        s0{i_sol}(:,k)     = mcoef{k}*Sjp';
        if Gapp == 1
            g0{i_sol}(:,k) = g0{i_sol}(:,k) - (mcoef{k}*Hfjp'*1e3)';
        end
    end
end
end