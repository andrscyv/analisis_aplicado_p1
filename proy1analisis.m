clear; clc; format long;
f = 'rosenbrock';
x0 = [3.5 4.5]';
[xf, iter, delta] = miregion(f, x0);