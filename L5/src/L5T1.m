clear
clc
close all

img_path = "img/";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end
%% Дано
syms t;

h = 2;
x = @(t) phi(t);
dxdt = @(t) -sign(x(t-h));
sim_time = 21.5;
%% Решение
t = linspace(0, sim_time, 1000);

p = figure;
plot(t, x_c(t));
grid on
xlim([0, sim_time])
xlabel('t')
ylabel('x')
print(p, img_path + "Z1" + ".png", '-dpng', '-r300');

p = figure;
plot(t, dxdt(t));
grid on
xlim([0, sim_time])
xlabel('t')
ylabel('$\dot{x}$', 'Interpreter', 'latex')
print(p, img_path + "Z2" + ".png", '-dpng', '-r300');


function p = phi(s)
    p = 0.5.*(s >= -2 & s < -1) + (-s - 0.5).*(s >= -1 & s <= 0)+...
        (-s - 0.5) .* (0   < s & s <= 1.5)+...
        ( s - 3.5) .* (1.5 < s & s <= 5.5  )+...
        (-s + 7.5) .* (5.5 < s & s <= 9.5  )+...
        ( s - 11.5) .* (9.5 < s & s <= 13.5  )+...
        (-s+ 15.5) .* (13.5 < s & s <= 17.5  )+...
        ( s - 19.5) .* (17.5 < s & s <= 21.5  );
end


function p = x_c(t)
    p = (-t - 0.5) .* (0   < t & t <= 1.5)+...
        ( t - 3.5) .* (1.5 < t & t <= 5.5  )+...
        (-t + 7.5) .* (5.5 < t & t <= 9.5  )+...
        ( t - 11.5) .* (9.5 < t & t <= 13.5  )+...
        (-t + 15.5) .* (13.5 < t & t <= 17.5  )+...
        ( t - 19.5) .* (17.5 < t & t <= 21.5  );
end