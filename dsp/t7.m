clc; clear; close all;

x = input('Enter sequence x[n] as a row vector: ');
h = input('Enter sequence h[n] as a row vector: ');

Lx = length(x);
Lh = length(h);
Ly = Lx + Lh - 1;
y_lin = zeros(1, Ly);

for n = 0:Ly-1
    kmin = max(0, n-(Lh-1));
    kmax = min(n, Lx-1);
    for k = kmin:kmax
        y_lin(n+1) = y_lin(n+1) + x(k+1) * h(n-k+1);
    end
end

disp('Linear Convolution Result y[n] = ');
disp(y_lin);