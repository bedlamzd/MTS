clear
clc
close all

rng(42);

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end

dxdt = @(x, r) [-x(1); r*x(2)- x(2)^3 + x(2)^5];

%% Якобиан
J = @(x, r) [
    -1, 0;
    0, r - 3*x(2)^2 + 5*x(2)^4
    ];

% Собственные числа
x = sym('x', [2 1]);
r = sym('r');

eig(J(x, r))

%% Положения равновесия

x2 = roots([1 0 -1 0 r 0]);
t = -1:0.001:1;
% for i=1:length(x2)
%     lambda = eig(J([0, x2(i)], r));
%     
%     figure;
%     hold on
%     grid on
%     t = -1:0.001:1;
%     q = eval(subs(lambda(2), r, t));
%     plot(t, real(q));
%     plot(t, imag(q));
%     plot(xlim, [lambda(1), lambda(1)]);
%     
%     
% %     r0 = solve(real(lambda(2)) == -1, r)
% %     if isempty(r0)
% %        continue 
% %     end
%     
% %     figure
% %     hold on
% %     phasePlot2(@(x) dxdt(x, 0.5), -1:0.1:1);
% %     for j = 1:10
% %         x0 = rand(2, 1)*(-1)^round(rand()*2);
% % %         for r_i = 0.25:0.1:1
% %             solPlot2(@(x) dxdt(x, 0.5), [0, 10], x0);
% % %         end
% %         xlim([-1 1]);
% %         ylim(xlim);
% %     end
% 
% end

x_ast = eval(subs(x2, r, t));

str = {};
for i = 1:length(x2)
    str{end+1} = ['$' latex(x2(i)) '$'];
end

figure
plot(t, real(x_ast));
grid on


legend(str, 'Interpreter', 'latex');

figure
plot(t, imag(x_ast));
grid on
legend(str, 'Interpreter', 'latex');


%% display params
width = 1;
height = 1;
n = 10;

%% (0, 0)
tspan = [0 10];

%% r = 0
r_i = 0;
center = [0, eval(subs(x2(1), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r0", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r0-lin", '-dpng', '-r300')

%% r = -1
r_i = -1;
center = [0, eval(subs(x2(1), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r-1", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r-1-lin", '-dpng', '-r300')

%% r = 1
r_i = 1;
center = [0, eval(subs(x2(1), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r1", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(1)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x1r1-lin", '-dpng', '-r300')

%% (0, x2_2) x45
%% r = 0
r_i = 0;
center = [0, eval(subs(x2(2), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(2)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x45r0", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(2)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x45r0-lin", '-dpng', '-r300')
%% r = 1/5 x45
r_i = 1/5;
center = [0, eval(subs(x2(2), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(2)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x45r02", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(2)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x45r02-lin", '-dpng', '-r300')


%% (0, x2_3) x23
%% r = 1/4
r_i = 1/4;
center = [0, eval(subs(x2(3), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(3)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x23r025", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(3)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x23r025-lin", '-dpng', '-r300')

%% r = 1/8 x23
r_i = 1/8;
center = [0, eval(subs(x2(3), r, r_i))];
[x_lims, y_lims] = getLimits(center, width, height);
lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
s = lims(1):diff(lims)/30:lims(2);
[x1_0, x2_0] = getInitials(center, width, height, n);

f = @(x) dxdt(x, r_i);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(3)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x23r0125", '-dpng', '-r300')

f = Linearize(@(x) J(x, r_i), center);

h = figure;

phasePlot2(f, s);
for i = 1:numel(x1_0)
    solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
end
ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(3)));
title(ttl,'Interpreter','latex');
sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
subtitle(sttl,'Interpreter','latex');
xlim(x_lims);
ylim(y_lims);

print(h, img_path + "x23r0125-lin", '-dpng', '-r300')


 %% r = -4 x23
 r_i = -4;
 center = [0, eval(subs(x2(3), r, r_i))];
 [x_lims, y_lims] = getLimits(center, width, height);
 lims = [min([x_lims, y_lims]), max([x_lims, y_lims])];
 s = lims(1):diff(lims)/30:lims(2);
 [x1_0, x2_0] = getInitials(center, width, height, n);

 f = @(x) dxdt(x, r_i);

 h = figure;

 phasePlot2(f, s);
 for i = 1:numel(x1_0)
     solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
 end
 ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$', r_i, center(2), latex(x2(3)));
 title(ttl,'Interpreter','latex');
 sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
 subtitle(sttl,'Interpreter','latex');
 xlim(x_lims);
 ylim(y_lims);

 print(h, img_path + "x23r-4", '-dpng', '-r300')

 f = Linearize(@(x) J(x, r_i), center);

 h = figure;

 phasePlot2(f, s);
 for i = 1:numel(x1_0)
     solPlot2(f, tspan, [x1_0(i), x2_0(i)]);
 end
 ttl = sprintf('$r = %0.2f; x^* = (0.00, %0.2f); x2 = %s$;', r_i, center(2), latex(x2(3)));
 title(ttl,'Interpreter','latex');
 sttl = sprintf('$\\lambda_1 = %0.2f; \\lambda_2 = %0.2f$', eig(J([0, center(2)], r_i)));
 subtitle(sttl,'Interpreter','latex');
 xlim(x_lims);
 ylim(y_lims);

 print(h, img_path + "x23r-4-lin", '-dpng', '-r300')

%% Utils
function phasePlot2(f, s)

    [x, y] = meshgrid(s, s);

    u = zeros(size(x));
    v = zeros(size(y));

    for i = 1:numel(x)
       dx = f([x(i), y(i)]);
       u(i) = dx(1);
       v(i) = dx(2);
    end
    
    h = gcf;
    hold on
    quiver(x, y, u, v, 2, 'b');
    xlabel("x1");
    ylabel("x2");
    axis equal;
    grid on

end

function solPlot2(f, tspan, x0)
    [~, y] = ode45(@(t, x) f(x), tspan, x0);

    plot(y(:,1), y(:,2), 'r');
    grid on
end

function [x_lims, y_lims] = getLimits(center, width, height)
    x_lims = [center(1) - width, center(1) + width];
    y_lims = [center(2) - height, center(2) + height];
end

function [x1, x2] = getInitials(center, width, height, n)
    [x_lims, y_lims] = getLimits(center, width, height);
    x1 = linspace(x_lims(1), x_lims(2), n);
    x2 = linspace(y_lims(1), y_lims(2), n);
    
    [x1, x2] = meshgrid(x1, x2);
end

function f_lin = Linearize(J, x_ast)
    f_lin = @(x) J(x_ast)*(x(:) - x_ast(:));
end
