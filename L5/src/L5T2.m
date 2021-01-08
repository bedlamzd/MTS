clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end

%% Дано
syms t x(t) tau(t) P q A A_1 eps;
assume(P > 0);
assume(q > 0);

psi = [
    A'*P + P*A + q*(1+eps)*P P*A_1;
    A_1'*P -q*P
    ];

latex(psi)

A = -2;
A_1 = -0.1;

eps = 1/2;
dx(t) = A*x(t) + A_1*x(t-tau(t));

psi = [
    A'*P + P*A + q*(1+eps)*P P*A_1;
    A_1'*P -q*P
    ];

latex(psi)

for n = 1:length(psi)
   minor = psi(1:n, 1:n);
   parity = (mod(n, 2) == 1);
   a = det(minor) * parity;
   b = det(minor) * ~parity;
   latex(det(minor))
   S = solve(a < b, q, 'ReturnConditions', true);
   vpa(S.conditions, 3)
end
