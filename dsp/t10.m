clc; clear; close all;

X = input('Enter DFT sequence X[k]: ');
N = length(X);
x_reconstructed = zeros(1, N);

for n = 0:N-1
    for k = 0:N-1
        x_reconstructed(n+1) = x_reconstructed(n+1) + X(k+1)*exp(1j*2*pi*k*n/N);
    end
    x_reconstructed(n+1) = x_reconstructed(n+1)/N;
end

x_reconstructed = real(x_reconstructed);

disp('Reconstructed sequence from IDFT:');
disp(x_reconstructed);

n = 0:N-1;
figure;
subplot(3,1,1);
stem(n, abs(X), 'filled'); title('Magnitude of input X[k]');
xlabel('k'); ylabel('|X[k]|');

subplot(3,1,2);
stem(n, angle(X), 'filled'); title('Phase of input X[k]');
xlabel('k'); ylabel('Angle');

subplot(3,1,3);
stem(n, x_reconstructed, 'filled'); title('IDFT Output x[n]');
xlabel('n'); ylabel('Amplitude');