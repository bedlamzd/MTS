rng(42);

clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
    mkdir(img_path)
end

%% Дана система

f = @(x) [
    -3*x(1)+x(2);...
    x(1)-4*x(2)-atan(2*x(2))
    ];

%% Линеаризация системы
A = @(x) [
    -3 1;
    1 -(6 + 16*x(2)^2)/(1 + 4*x(2)^2)
];

x0 = rand(2, 1); % Начальное приближение для решения f(x) = 0
x_ast = fsolve(f, x0); % найдём одно из положений равновесия

f_lin = @(x) A(x_ast) * (x - x_ast); % линеаризованая система


%% Исследование положения равновесия
lambda = eig(A(x_ast));

disp(lambda); % оба корня действительные и отрицательные, 
              % что значит равновесие типа "узел"

s = linspace(-10, 10, 15);

[x, y] = meshgrid(s, s);

u = zeros(size(x));
v = zeros(size(y));

for i = 1:numel(x)
   dx = f([x(i), y(i)]);
   u(i) = dx(1);
   v(i) = dx(2);
end

h = figure;
quiver(x, y, u, v, 3, 'b');
ax = h.CurrentAxes;
set(ax, 'XAxisLocation', 'origin');
set(ax, 'YAxisLocation', 'origin');
title("Фазовая диаграмма системы");
xlabel("x1");
ylabel("x2");
axis tight equal;
grid on

print(h, img_path + "phase", '-dpng', '-r300')

%% Построим графики систем
t0 = 0;
tf = 10;
tspan = [t0 tf]; % пределы графика

%% Исходная система
[t, y] = ode45(@(t, x) f(x), tspan, x0);

h = figure;
plot(t, y);
title("Исходная система");
subtitle(sprintf('x_0 = [%0.3f %0.3f]', x0));
xlabel('t');
ylabel('x');
legend('x_1', 'x_2');
grid on

print(h, img_path + "initial", '-dpng', '-r300')

%% Линеаризованая система
[t, y] = ode45(@(t, x) f_lin(x), tspan, x0);

h = figure;
plot(t, y);
title("Линеаризованая система");
subtitle(sprintf('x_0 = [%0.3f %0.3f]', x0));
xlabel('t');
ylabel('x');
legend('x_1', 'x_2');
grid on

print(h, img_path + "linearized", '-dpng', '-r300')
