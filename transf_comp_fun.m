function X = transf_comp_fun(x_old,new_comp2old_comp)
% Component transformation
% The matrix of compositions of the old components and new components
A = new_comp2old_comp;
% The vector of composition in terms of old component
x = x_old';
% Check number of independent components
neq = rank(A);
% Make the set of equations, by using the rank
ieq = 0;
B = A(1,:);
b = x(1,:);
for i = 2:size(x,1) % Keep adding equations until rank is reached
    if rank([B; A(i,:)]) > rank(B)
        ieq = ieq + 1;
        B = [B; A(i,:)];
        b = [b; x(i,:)];
    end
    if rank(B) == neq, break,end
end
% Solve system of equations
X = B\b;
% Remove numerical noise
X( X> -1e-8 & X<  1e-8 ) = 0;
X( X>1-1e-8 & X<1+1e-8 ) = 1;