clear
clc
close all

img_path = ".\..\img\";

if ~exist(img_path, 'dir')
    mkdir(img_path)
end

A = [
    0 1 0;
    0 0 1;
    0.1 0.3 0.2
    ];
b = [
    0.1;
    -0.3;
    0.2
    ];
C = [0 0 1];

syms s;

W = C*(s * eye(3) - A)^-1*b;

