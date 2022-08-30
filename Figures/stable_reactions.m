function alph = stable_reactions(A_react,dG,solid_names)
rxn = zeros(size(dG,1)*size(dG,2),size(A_react,1));
ass = zeros(size(dG,1)*size(dG,2),size(A_react,1),size(A_react,2));
alph = zeros(size(dG,1)*size(dG,2),length(solid_names));
% Make phase alpha matrix
for ir = 1:size(A_react,1)
    rxn(dG(:,:,ir)<0,ir) = -1;
    rxn(dG(:,:,ir)>0,ir) = 1;
    ass(rxn(:,ir)==-1,ir,A_react(ir,:)>0) = 1;
    ass(rxn(:,ir)== 1,ir,A_react(ir,:)<0) = 1;
end
for ip = 1:length(solid_names)
    irxn = A_react(:,ip)~=0; % index of reactions in which the phase occurs
    nrxn = sum(irxn); % number of reactions the phase occurs
    alph(sum(ass(:,irxn,ip),2)==nrxn,ip) = 1; % if it occurs in all reactions it is not excluded (I hope it means it is stable)0
end