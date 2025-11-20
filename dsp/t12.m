clc; clear all; close all;

function X = ifft2pt(x)
    X = [x(1)+x(2), x(1)-x(2)]/2;
end

function X = ifft4pt(x)
    G = ifft2pt([x(1) x(3)]);
    H = ifft2pt([x(2) x(4)]);
    X = [G(1)+H(1), G(2)+1j*H(2), G(1)-H(1), G(2)-1j*H(2)]/2;
end

function X = ifft8pt(x)
    G = ifft4pt([x(1) x(3) x(5) x(7)]);
    H = ifft4pt([x(2) x(4) x(6) x(8)]);
    W = exp(1j*2*pi*(0:3)/8);
    X = [G+H.*W, G-H.*W]/2;
end

xk = input('Enter input: ');
n = length(xk);

if n==2
    xn = ifft2pt(xk);
elseif n==4
    xn = ifft4pt(xk);
elseif n==8
    xn = ifft8pt(xk);
else
    error('Only 2/4/8 supported');
end

disp('IFFT:');
disp(xn);

subplot(3,1,1), stem(abs(xk)), title('|X(k)|')
subplot(3,1,2), stem(angle(xk)), title('Phase')
subplot(3,1,3), stem(real(xn)), title('Output x(n)')