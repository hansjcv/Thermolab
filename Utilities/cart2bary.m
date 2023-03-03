function [x2d,y2d] = cart2bary(X1_2d,X2_2d)
a    = 2/3*sqrt(3);%1;
b    = 1/3*sqrt(3);%1/2;
A    = [a b;
    0 1];
xc = A*[X1_2d(:) X2_2d(:)]';
x = xc(1,:);
y = xc(2,:);
x2d = reshape(x,size(X1_2d,1),size(X1_2d,2));
y2d = reshape(y,size(X1_2d,1),size(X1_2d,2));