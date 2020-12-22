clear
clc
close all

A = [-2 1;...
     1 -2];
B = [0;
     1];
C = [1;
     0];


sigma = @(x) c*x;
ksi = @(sigma) (exp(sigma) - exp(-sigma))/(exp(sigma) + exp(-sigma));
dxdt = @(x) A*x + B*ksi;

mu1 = 0;
mu2 = 1;

lambda = eig(A);

disp(lambda);

% Выбирая mu0 = mu1 = 0 система, при ksi = mu0*sigma, принимает вид 
% dxdt = A*x, а т.к. собственные числа матрицы действительны и < 0, 
% то система устойчива

syms s
w = sym('w', 'real');
W(s) = C'*(s*eye(length(A)) - A)^-1*B;
W_0(s) = -W(s);

real((1+mu1*W_0(1i*w))*conj(1+mu2*W_0(1i*w)))