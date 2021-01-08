clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end
%% Дано
syms t;

h = 2;
% x = @(t) phi(t);
dxdt = @(t) -sign(x(t-h));

%% Решение
t = linspace(-h, 5*h, 1000);

p = figure;
plot(t, x(t));
grid on
xlim([-2, 6])
xlabel('t')
ylabel('x')
print(p, img_path + "Z1" + ".png", '-dpng', '-r300');

p = figure;
plot(t, dxdt(t));
grid on
xlim([-2, 6])
xlabel('t')
ylabel('$\dot{x}$', 'Interpreter', 'latex')

function p = phi(t)
    p = 0.5.*(t >= -2 & t < -1) + (-t - 0.5).*(t >= -1 & t <= 0);
end

function p = x(t)
    p = phi(t) +...
        (-t - 0.5) .* (0   < t & t <= 1.5)+...
        ( t - 3.5) .* (1.5 < t & t <= 2  )+...
        ( t - 3.5) .* (2   < t & t <= 4  )+...
        ( t - 3.5) .* (4   < t & t <= 5.5)+...
        (-t + 7.5) .* (5.5 < t & t <= 6);
end