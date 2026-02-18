clear,clf,addpath ../ ../EOS
% Example phase diagrams Fig 4
T = linspace(400,850,100) + 273.15;
P = linspace(0.1,6,101)*1e9;            
X    = {'Si','Mg','H','O'};
load tl_dataset
% phases = {'atg,tc-ds55','br,tc-ds55','fo,tc-ds55','ta,tc-ds55','en,tc-ds55','anth,tc-ds55','H2O,tc-ds55'};
phases = {'atg,tc-ds633','br,tc-ds633','fo,tc-ds633','ta,tc-ds633','en,tc-ds633','anth,tc-ds633','H2O,tc-ds633'};
[T2d,P2d] = ndgrid(T,P);
[G,Nphs] = tl_gibbs_energy(T2d(:),P2d(:),phases);
Nsys = Nphs(:,1); 
LB  = zeros(1,size(G,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iPT = 1:length(T2d(:))
    alph(iPT,:) = linprog(G(:,iPT),[],[],Nphs,Nsys,LB); % The Gibbs energy minimization       
end
% Plotting
assemblage_id = zeros(length(T)*length(P),length(phases));
for i = 1:size(assemblage_id,1)
    assemblage_id(i,1:length(find(alph(i,:)>0))) = find(alph(i,:)>0);
end
% Plot phase diagram section
figure(1)
tl_psection(T-273.15,P/1e9,X,assemblage_id,phases);
xlabel('T (\circC)'),ylabel('P(GPa)')
set(gca,'FontSize',12)