clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
   mkdir(img_path);
end

%% Дано
initial_nom = [1 2];
initial_den = [1 4 5];

syms s h z;

W = makeW(initial_nom, initial_den);
W_eul = subs(W, s, (s-1)/h);
W_tas = subs(W, s, 2/h*(s-1)/(s+1));

prettyLatex(subs(W_eul, s, z), z)
prettyLatex(subs(W_tas, s, z), z)

A = poly2sym(initial_den, s);
h_est = eval(min(2*abs(real(solve(A)))./abs(solve(A)).^2));

p = figure;
stepplot(sym2tfs(W));
title('W');
grid on
print(p, img_path + "init.png", '-dpng', '-r300');

for h_i = [h_est/2, h_est, h_est*2, h_est*100]
    p = figure;
    stepplot(sym2tfz(subs(W_eul, h, h_i), h_i));
    grid on
    title(sprintf('W_{eul} h = %0.2f', h_i));
    print(p, img_path + "eul" + sprintf('h%g', h_i) + ".png",...
     '-dpng', '-r300');
    
    p = figure;
    stepplot(sym2tfs(W));
    hold on
    stepplot(sym2tfz(subs(W_eul, h, h_i), h_i));
    grid on
    title(sprintf('W_{eul} h = %0.2f', h_i));
    legend("Исходная", "Дискретизация по Эйлеру");
    print(p, img_path + "init-eul" + sprintf('h%g', h_i) + ".png",...
     '-dpng', '-r300');
    
end

for h_i = [h_est/2, h_est, h_est*2, h_est*100]
    p = figure;
    stepplot(sym2tfz(subs(W_tas, h, h_i), h_i));
    grid on
    title(sprintf('W_{tas}; h = %0.2f', h_i));
    print(p, img_path + "tas" + sprintf('h%g', h_i) + ".png",...
     '-dpng', '-r300');
    
    p = figure;
    stepplot(sym2tfs(W));
    hold on
    stepplot(sym2tfz(subs(W_tas, h, h_i), h_i));
    grid on
    title(sprintf('W_{tas}; h = %0.2f', h_i));
    legend("Исходная", "Дискретизация по Тастину");
    print(p, img_path + "init-tas" + sprintf('h%g', h_i) + ".png",...
     '-dpng', '-r300');
end

function W = makeW(nom, den)
    syms s;
    W = poly2sym(nom, s)/poly2sym(den, s);
end

function TF = sym2tfs(W)
    [num, den] = numden(W);
    tfn = sym2poly(num);
    tfd = sym2poly(den);
    TF = tf(tfn, tfd);
end

function TF = sym2tfz(W, delay)
    [num, den] = numden(W);
    tfn = sym2poly(num);
    tfd = sym2poly(den);
    TF = tf(tfn, tfd, delay);
end

function str = prettyLatex(W, s)
    [num, den] = numden(W);
    num = simplify(collect(num, s));
    den = simplify(collect(den,s));
    str = latex(num/den);
end
