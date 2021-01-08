clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end

%% Дано
A = [
    3 1;
    2 2
    ];
A_1 = [
    -5 1;
    -2 -4
    ];

dx = @(t, x, Z) A * x + A_1*Z;
phi = @(t) [5*sin(t); 2*cos(t)+1];

%% Решение
MaxTime = 15;
MinTime = 0;
tspan = [0 MaxTime];
tgrid = linspace(MinTime, MaxTime, 1000);

for h = [0.15 0.1770 0.2]
    sol = dde23(dx, h, phi, tspan);
    x = deval(sol, tgrid);

    p = figure;
    plot(tgrid, x);
    grid on
    title(sprintf('h = %.2f', h));
    legend('x_1','x_2');
    xlabel('t');
    ylabel('x(t)');
    print(p, img_path + "L6" + sprintf('h%.2f', h) + ".png", '-dpng', '-r300');
    
end

h = 0.1770;
[~, P, R] = solveDesc(A, A_1, h);
[~, K] = solveLap(A, A_1, value(P));

value(P)
value(R)
value(K)

MaxTime = 2;
MinTime = 0;
tspan = [0 MaxTime];
tgrid = linspace(MinTime, MaxTime, 1000);
dx = @(t, x, Z) (A+K)*x + A_1*Z;
for h = [0.15 0.1770 0.2]
    sol = dde23(dx, h, phi, tspan);
    x = deval(sol, tgrid);

    p = figure;
    plot(tgrid, x);
    grid on
    title(sprintf('h = %.2f', h));
    legend('x_1','x_2');
    xlabel('t');
    ylabel('x(t)');
    print(p, img_path + "L62" + sprintf('h%.2f', h) + ".png", '-dpng', '-r300');
    
end

function [sol, K] = solveLap(A, A1, P)
    AplusK = sdpvar(2,2);
    Q = sdpvar(2,2);

    TH=blkvar;

    TH(1,1) = ((AplusK)')*P + P*AplusK+Q;
    TH(1,2) = (A1')*P;
    TH(2,2) = -Q;
    TH=sdpvar(TH);

    F = [TH<=-0.0001, Q>=0.0001];

    options=sdpsettings('solver','sedumi','verbose',0);

    sol=optimize(F)

    AplusK = value(AplusK);
    K = AplusK - A;
end

function [sol, P, R] = solveDesc(A, A1, h)
    n=length(A);        

    Q = A+A1;

    P = sdpvar(n);
    R = sdpvar(n);
    P2 = sdpvar(n);
    P3 = sdpvar(n);

    TH = blkvar;

    TH(1,1) = (P2')*Q+(Q')*P2;
    TH(1,2) = P - P2'+(Q')*P3;
    TH(1,3) = -h*(P2')*A1;
    TH(2,2) = -P3-P3'+h*R;
    TH(2,3) = -h*(P3')*A1;
    TH(3,3) = -h*R;

    TH=sdpvar(TH);

    F=[TH<=-1e-4, R>=1e-4, P>=1e-4];

    options=sdpsettings('solver','sedumi','verbose',0);

    sol=optimize(F)

end