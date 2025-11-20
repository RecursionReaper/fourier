clc; clear; close all;

x = input('Enter sequence x[n] as a row vector: ');
h = input('Enter sequence h[n] as a row vector: ');
N = input('Enter N for N-point circular convolution: ');

x_pad = x;
h_pad = h;

x_pad(end+1:N) = 0;
h_pad(end+1:N) = 0;

x_pad = x_pad(1:N);
h_pad = h_pad(1:N);

y_circ = zeros(1, N);

for n = 0:N-1
    for k = 0:N-1
        y_circ(n+1) = y_circ(n+1) + x_pad(k+1) * h_pad(mod(n-k, N)+1);
    end
end

disp(['Circular Convolution Result y[n] for N = ' num2str(N) ':']);
disp(y_circ);