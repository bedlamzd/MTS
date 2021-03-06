rng(42); % репродуцируемые случайные значения

clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
    mkdir(img_path)
end

%% Дано
A = [-3 1 5;...
     3 -4 1;...
     3 0 -4];
B = eye(3);

%% Моделирование 
x0 = rand(3, 1); % ненулевые начальные условия
t0 = 0;
tf = 10;
tspan = [t0 tf]; % время моделирования

%% С ненулевыми начальными условиями без управления
dxdt = @(t, x) A*x; % моделируемая система

[t, y] = ode45(dxdt, tspan, x0);

h = figure;
plot(t, y);
title("Без управления");
subtitle(sprintf('x_0 = [%0.3f %0.3f %0.3f]', x0));
xlabel('t');
ylabel('x');
legend('x_1', 'x_2', 'x_3');
grid on

print(h, img_path + "no_control", '-dpng', '-r300')

lambda_open = eig(A);

%% Те же начальные условия, но с управлением 
% (ассимптотическая устойчивость)
k = -max(real(eig(A))); % граничный коэффициент усиления k*
u = @(x) k*x; % управление системы
dxdt = @(t, x) A*x + B*u(x); % моделируемая система

[t, y] = ode45(dxdt, tspan, x0);

h = figure;
plot(t, y);
title("С управлением (Ассимптотическая устойчивость)");
subtitle(sprintf('x_0 = [%0.3f %0.3f %0.3f]; k=%0.3f', x0, k));
xlabel('t');
ylabel('x');
legend('x_1', 'x_2', 'x_3');
grid on

print(h, img_path + "assymp_control", '-dpng', '-r300')

lambda_closed = eig(A + k*B);
