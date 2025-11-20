clc; clear all; close all;

function X = fft2pt(x)
    X = [x(1)+x(2), x(1)-x(2)];
end

function X = fft4pt(x)
    G = fft2pt([x(1) x(3)]);
    H = fft2pt([x(2) x(4)]);
    X = [G(1)+H(1), G(2)-1j*H(2), G(1)-H(1), G(2)+1j*H(2)];
end

function X = fft8pt(x)
    G = fft4pt([x(1) x(3) x(5) x(7)]);
    H = fft4pt([x(2) x(4) x(6) x(8)]);
    W = exp(-1j*2*pi*(0:3)/8);
    X = [G+H.*W, G-H.*W];
end

x = input('Enter input: ');
n = length(x);

if n==2
    X = fft2pt(x);
elseif n==4
    X = fft4pt(x);
elseif n==8
    X = fft8pt(x);
else
    error('Only 2/4/8 supported');
end

disp('FFT:');
disp(X);

subplot(3,1,1), stem(real(x)), title('Input x(n)')
subplot(3,1,2), stem(abs(X)), title('|X(k)|')
subplot(3,1,3), stem(angle(X)), title('Phase')