clear
clc
close all

A = [-1 -3;...
     -1 -4];
B = [1;
     0];
C = [0;
     1];


sigma = @(x) c*x;
ksi = @(sigma) (abs(sigma) < 1)*2*sigma + (abs(sigma) >= 1)*2*sign(sigma);
dxdt = @(x) A*x + B*ksi;

mu0 = 2;

lambda = eig(A);

disp(lambda);

syms s
w = sym('w', 'real');
nu = sym('nu', 'real');
W(s) = C'*(s*eye(length(A)) - A)^-1*B;
W_0(s) = -W(s);

mu0^(-1) + real((1+1i*w*nu)*W_0(1i*w))