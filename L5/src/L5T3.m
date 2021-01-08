clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end

%% Дано
A = [
    -4 1;
    -2 -4
    ];
A_1 = [
    1 2;
    1 -1
    ];
h = 3;

dx = @(t, x, Z) A * x + A_1*Z;
phi = @(t) [5*sin(t); 2*cos(t)+1];

MaxTime = 15;
MinTime = 0;
tspan = [0 MaxTime];
tgrid = linspace(MinTime, MaxTime, 1000);

% Solving
control = false;
sol = dde23(dx, h, phi, tspan);
x = deval(sol, tgrid);

p = figure;
plot(tgrid, x);
grid on
legend('x_1','x_2');
xlabel('t');
ylabel('x(t)');
print(p, img_path + "Z3" + ".png", '-dpng', '-r300');

[P, Q] = solveIneq(A, A_1);

value(P)
value(Q)

function [P, Q] = solveIneq(A, A1)
    n=length(A);            

    P=sdpvar(n);   
    Q=sdpvar(n);
    TH=blkvar;       
    TH(1,1)=A'*P+P*A+Q;
    TH(1,2)=P*A1;
    TH(2,2)=-Q;
    TH=sdpvar(TH);

    F=[TH<=-0.0001,P>=0.0001,Q>=0.0001];
    options=sdpsettings('solver','sedumi','verbose',0);
    sol=optimize(F)
end
