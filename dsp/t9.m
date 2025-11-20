clc; clear; close all;

x = input('Enter input sequence x[n]: ');
N = length(x);
X = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1)*exp(-1j*2*pi*k*n/N);
    end
end

disp('DFT of x:');
disp(X);

k = 0:N-1;
figure;
subplot(3,1,1);
stem(0:N-1, x, 'filled'); title('Input Sequence x[n] for DFT');
xlabel('n'); ylabel('Amplitude');

subplot(3,1,2);
stem(k, abs(X), 'filled'); title('Magnitude Spectrum of x[k]');
xlabel('k'); ylabel('|X[k]|');

subplot(3,1,3);
stem(k, angle(X), 'filled'); title('Phase Spectrum of x[k]');
xlabel('k'); ylabel('LX[k]');