[20/11/25, 2:57:36 AM] Cocomelon: 1) Generation of standard CT signal






[20/11/25, 2:57:40 AM] Cocomelon:



 t = -10:0.1:10;
impulse_c = (t==0);
step_c    = (t>=0);
ramp_c    = t.*(t>=0);
expo_c    = exp(t);

figure;
subplot(2,2,1); plot(t,impulse_c); title('Continuous Impulse δ(t)'); xlabel('t'); 
ylabel('Amplitude');
subplot(2,2,2); plot(t,step_c);    title('Continuous Step u(t)');    xlabel('t'); 
ylabel('Amplitude');
subplot(2,2,3); plot(t,ramp_c);    title('Continuous Ramp r(t)');    xlabel('t'); 
ylabel('Amplitude');
subplot(2,2,4); plot(t,expo_c);    title('Continuous Exponential e^t'); 
xlabel('t'); ylabel('Amplitude');
sgtitle('Figure 3: Continuous Elementary Signals');

t = -10:0.1:10; w = pi/3;
yc1 = sin(w*t);
yc2 = cos(w*t);
yc3 = sin(w*t)./(w*t);
yc3(t==0) = 1;

figure;
subplot(3,1,1); plot(t,yc1); 
title('Continuous Sine sin(\omega t)'); xlabel('t'); ylabel('Amplitude');

subplot(3,1,2); plot(t,yc2); 
title('Continuous Cosine cos(\omega t)'); xlabel('t'); ylabel('Amplitude');

subplot(3,1,3); plot(t,yc3); 
title('Continuous Sinc sinc(\omega t)'); xlabel('t'); ylabel('Amplitude');

sgtitle('Figure 4: Continuous Trigonometric Signals');
[20/11/25, 2:58:11 AM] Cocomelon:



 2) Generation of standard DT signal
[20/11/25, 2:58:13 AM] Cocomelon:



 n = -10:10;
impulse = (n==0);
step    = (n>=0);
ramp    = n.*(n>=0);
expo    = exp(n);

figure;
subplot(2,2,1); stem(n,impulse); title('Discrete Impulse δ[n]'); xlabel('n'); 
ylabel('Amplitude');
subplot(2,2,2); stem(n,step);    title('Discrete Step u[n]');    xlabel('n'); 
ylabel('Amplitude');
subplot(2,2,3); stem(n,ramp);    title('Discrete Ramp r[n]');    xlabel('n'); 
ylabel('Amplitude');
subplot(2,2,4); stem(n,expo);    title('Discrete Exponential e^n'); xlabel('n'); 
ylabel('Amplitude');
sgtitle('Figure 1: Discrete Elementary Signals');

n = -10:10; w = pi/3;
y1 = sin(w*n);
y2 = cos(w*n);
y3 = sin(w*n)./(w*n);
y3(n==0) = 1;

figure;
subplot(3,1,1); stem(n,y1,'filled'); 
title('Discrete Sine sin(\omega n)'); xlabel('n'); ylabel('Amplitude');

subplot(3,1,2); stem(n,y2,'filled'); 
title('Discrete Cosine cos(\omega n)'); xlabel('n'); ylabel('Amplitude');

subplot(3,1,3); stem(n,y3,'filled'); 
title('Discrete Sinc sinc(\omega n)'); xlabel('n'); ylabel('Amplitude');

sgtitle('Figure 2: Discrete Trigonometric Signals');





[20/11/25, 3:12:41 AM] Cocomelon: 3) Generation of customised DT and CT signal


clc; clear; close all;

x = input('Enter sequence x: ')

n = input('Enter index vector n: ')

a = input('Enter scaling factor a: ')

m = input('Enter shift value m: ')

figure(1);
subplot(2,2,1); stem(n, x, 'filled'); title('x[n]'); xlabel('n'); 
ylabel('Amplitude'); grid on;
subplot(2,2,2); stem(n + m, x, 'filled'); title('x[n - m] (Shift right by +m)'); 
xlabel('n'); ylabel('Amplitude'); grid on;
subplot(2,2,3); stem(-n, x, 'filled'); title('x[-n] (Folded)'); xlabel('n'); 
ylabel('Amplitude'); grid on;
subplot(2,2,4); stem(n / a, x, 'filled'); title('x[a n] (Time-scaling)'); 
xlabel('n'); ylabel('Amplitude'); grid on;




[20/11/25, 3:12:52 AM] Cocomelon: 4) an+b
[20/11/25, 3:13:51 AM] Cocomelon: 




clc; clear; close all;

x = input('Enter sequence x: ')

n = input('Enter index vector n: ')

a = input('Enter scaling factor a: ')

m = input('Enter shift value m: ')



n1 = n - m;
figure(1);
stem(n1 / a, x .* (mod(n1, a) == 0), 'filled'); 
‎Read more
[20/11/25, 3:14:13 AM] Cocomelon: 5) plot even part of DT signal
[20/11/25, 3:16:43 AM] Cocomelon: clc; clear; close all;

x = input('Enter the sequence x[n]: ');
n = input('Enter the corresponding indices n: ');

x_folded = fliplr(x);
n_folded = -fliplr(n);

xe = (x + x_folded) / 2;

disp('Original sequence x[n]:');
disp(x);

disp('Folded sequence x[-n]:');
disp(x_folded);

disp('Even part xe[n]:');
disp(xe);

figure;
subplot(3,1,1);
stem(n, x, 'filled'); grid on;
xlabel('n'); ylabel('x[n]');
title('Original Signal x[n]');

subplot(3,1,2);
stem(n_folded, x_folded, 'filled'); grid on;
xlabel('n'); ylabel('x[-n]');
title('Folded Signal x[-n]');

subplot(3,1,3);
stem(n, xe, 'filled'); grid on;
xlabel('n'); ylabel('xe[n]');
title('Even Part xe[n]');




[20/11/25, 3:16:55 AM] Cocomelon: 6) odd part of DT signal
[20/11/25, 3:16:56 AM] Cocomelon: 



clc; clear; close all;

x = input('Enter the sequence x[n]: ');
n = input('Enter the corresponding indices n: ');

x_folded = fliplr(x);
n_folded = -fliplr(n);

xo = (x - x_folded) / 2;

disp('Original sequence x[n]:');
disp(x);

disp('Folded sequence x[-n]:');
disp(x_folded);

disp('Odd part xo[n]:');
disp(xo);

figure;
subplot(3,1,1);
stem(n, x, 'filled'); grid on;
xlabel('n'); ylabel('x[n]');
title('Original Signal x[n]');

subplot(3,1,2);
stem(n_folded, x_folded, 'filled'); grid on;
xlabel('n'); ylabel('x[-n]');
title('Folded Signal x[-n]');

subplot(3,1,3);
stem(n, xo, 'filled'); grid on;
xlabel('n'); ylabel('xo[n]');
title('Odd Part xo[n]');



[20/11/25, 3:19:50 AM] Cocomelon: 7) linear convolution
[20/11/25, 3:19:51 AM] Cocomelon:



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





[20/11/25, 3:20:02 AM] Cocomelon: 8) circular convolution
[20/11/25, 3:20:03 AM] Cocomelon: 


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



[20/11/25, 3:28:11 AM] Cocomelon: 9) DFT of given seque
[20/11/25, 3:28:11 AM] Cocomelon:

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
xlabel('k'); ylabel('∠X[k]');


[20/11/25, 3:28:26 AM] Cocomelon: 10) IDFT of given seque
[20/11/25, 3:28:27 AM] Cocomelon: 



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





[20/11/25, 3:29:50 AM] Cocomelon: 11) DIT FFT
[20/11/25, 3:29:51 AM] Cocomelon: 



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




[20/11/25, 3:29:57 AM] Cocomelon: 12) DIT IFFT
[20/11/25, 3:30:02 AM] Cocomelon:




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




[20/11/25, 3:33:17 AM] Cocomelon: 13) Obtain frequency response of H(z) and pole zero p
[20/11/25, 3:33:22 AM] Cocomelon: 



num = input("Enter nr coeff: ");
den = input("Enter dr coeff: ");
[H,w]=freqz(num,den,512);
subplot(3,1,1);
plot(w,abs(H)); grid on;
title('Mag response');
subplot(3,1,2);
plot(w, angle(H)); grid on;
title('Phase response');
subplot(3,1,3);
zplane(num,den);
title('Pole zero plot');




[20/11/25, 3:33:30 AM] Cocomelon: 14) Obtain frequency response of h(n)and pole zero p
[20/11/25, 3:33:34 AM] Cocomelon:



 h = input("Enter the sequence h(n): ");
[H,w] = freqz(h, 1, 512);
subplot(3,1,1);
plot(w, abs(H)); grid on;
title('Magnitude Response of h(n)');
subplot(3,1,2);
plot(w, angle(H)); grid on;
title('Phase Response of h(n)');
subplot(3,1,3);
zplane(h, 1);
title('Pole-Zero Plot of h(n)');




[20/11/25, 3:37:34 AM] Cocomelon: 15) frequency samp
[20/11/25, 3:37:36 AM] Cocomelon: band pass
[20/11/25, 3:37:37 AM] Cocomelon:


 N = input('N: ');
wc1 = input('wc1: ');
wc2 = input('wc2: ');

fprintf('Band Pass Filter -> N = %d, wc1 = %.2f, wc2 = %.2f\n', N, wc1, wc2);

k = 0:N-1;
wk = 2*pi*k/N;

Hk = zeros(1, N);
Hk((wk >= wc1 & wk <= wc2) | (wk >= (2*pi - wc2) & wk <= (2*pi - wc1))) = 1;

h = real(ifft(Hk));
disp(h);

figure;
subplot(2,1,1);
stem(k, abs(Hk));
title('H(k) BPF'); grid on

[H, w] = freqz(h, 1, 512);
subplot(2,1,2);
plot(w, abs(H));
title('|H(\omega)| BPF'); grid on




[20/11/25, 3:38:01 AM] Cocomelon: low p
[20/11/25, 3:38:01 AM] Cocomelon: N = input('N: ');
wc1 = input('wc1: ');
wc2 = input('wc2: ');

fprintf('Low Pass Filter -> N = %d, wc1 = %.2f, wc2 = %.2f\n', N, wc1, wc2);

k = 0:N-1;
wk = 2*pi*k/N;

Hk = zeros(1, N);
Hk((wk >= wc1 & wk <= wc2) | (wk >= (2*pi - wc2) & wk <= (2*pi - wc1))) = 1;

h = real(ifft(Hk));
disp(h);

figure;
subplot(2,1,1);
stem(k, abs(Hk));
title('H(k) LPF'); grid on

[H, w] = freqz(h, 1, 512);
subplot(2,1,2);
plot(w, abs(H));
title('|H(\omega)| LPF'); grid on
[20/11/25, 3:38:09 AM] Cocomelon: high p
[20/11/25, 3:38:10 AM] Cocomelon: N = input('N: ');
wc1 = input('wc1: ');
wc2 = input('wc2: ');

fprintf('High Pass Filter -> N = %d, wc1 = %.2f, wc2 = %.2f\n', N, wc1, wc2);

k = 0:N-1;
wk = 2*pi*k/N;

Hk = zeros(1, N);
Hk((wk >= wc1 & wk <= wc2) | (wk >= (2*pi - wc2) & wk <= (2*pi - wc1))) = 1;

h = real(ifft(Hk));
disp(h);

figure;
subplot(2,1,1);
stem(k, abs(Hk));
title('H(k) HPF'); grid on

[H, w] = freqz(h, 1, 512);
subplot(2,1,2);
plot(w, abs(H));
title('|H(\omega)| HPF'); grid on
[20/11/25, 3:41:29 AM] Cocomelon: band 
[20/11/25, 3:41:29 AM] Cocomelon: N = input('N: ');
wc1 = input('wc1: ');
wc2 = input('wc2: ');

fprintf('Band Stop Filter -> N = %d, wc1 = %.2f, wc2 = %.2f\n', N, wc1, wc2);

k = 0:N-1;
wk = 2*pi*k/N;

Hk = zeros(1, N);
Hk(~((wk >= wc1 & wk <= wc2) | (wk >= (2*pi - wc2) & wk <= (2*pi - wc1)))) = 1;

h = real(ifft(Hk));
disp(h);

figure;
subplot(2,1,1);
stem(k, abs(Hk));
title('H(k) BSF'); grid on

[H, w] = freqz(h, 1, 512);
subplot(2,1,2);
plot(w, abs(H));
title('|H(\omega)| BSF'); grid on





[20/11/25, 3:41:41 AM] Cocomelon: 16) window
[20/11/25, 3:42:32 AM] Cocomelon: low p
[20/11/25, 3:42:32 AM] Cocomelon: clc; clear; close all
N = input('Enter the number of elements: ')
wc1 = input('Enter the lower cutoff frequency: ')
wc2 = input('Enter the higher cutoff frequency: ')

a = (N - 1) / 2;
n = -a:a;
mid = a + 1;

h = (sin(wc2*n) - sin(wc1*n)) ./ (pi*n);
h(mid) = (wc2 - wc1) / pi;

disp('Impulse Response h(n):');
disp(h);

h_r  = h .* 1
h_hn = h .* (0.5 + 0.5 * cos((pi*n) ./ a))
h_hm = h .* (0.54 + 0.46 * cos((pi*n) ./ a))
h_bl = h .* (0.42 + 0.5 * cos((pi*n) ./ a) + 0.08 * cos((2*pi*n) ./ a))

[Hr, w] = freqz(h_r, 1, 1024);
[Hn, ~] = freqz(h_hn, 1, 1024);
[Hm, ~] = freqz(h_hm, 1, 1024);
[Hb, ~] = freqz(h_bl, 1, 1024);

figure('Name','Band Pass FIR Filter');
subplot(2,2,1); plot(w, abs(Hr)); title('Rectangular Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,2); plot(w, abs(Hn)); title('Hanning Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,3); plot(w, abs(Hm)); title('Hamming Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,4); plot(w, abs(Hb)); title('Blackman Window'); xlabel('w'); ylabel('h(n)'); grid("on");


[20/11/25, 3:42:40 AM] Cocomelon: high p
[20/11/25, 3:42:40 AM] Cocomelon: 




clc; clear; close all
N = input('Enter the number of elements: ')
wc1 = input('Enter the lower cutoff frequency: ')
wc2 = input('Enter the higher cutoff frequency: ')

a = (N - 1) / 2;
n = -a:a;
mid = a + 1;

h = (sin(wc2*n) - sin(wc1*n)) ./ (pi*n);
h(mid) = (wc2 - wc1) / pi;

disp('Impulse Response h(n):');
disp(h);

h_r  = h .* 1;
h_hn = h .* (0.5 + 0.5 * cos((pi*n) ./ a))
h_hm = h .* (0.54 + 0.46 * cos((pi*n) ./ a))
h_bl = h .* (0.42 + 0.5 * cos((pi*n) ./ a) + 0.08 * cos((2*pi*n) ./ a))

[Hr, w] = freqz(h_r, 1, 1024);
[Hn, ~] = freqz(h_hn, 1, 1024);
[Hm, ~] = freqz(h_hm, 1, 1024);
[Hb, ~] = freqz(h_bl, 1, 1024);

figure('Name','Band Pass FIR Filter');
subplot(2,2,1); plot(w, abs(Hr)); title('Rectangular Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,2); plot(w, abs(Hn)); title('Hanning Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,3); plot(w, abs(Hm)); title('Hamming Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,4); plot(w, abs(Hb)); title('Blackman Window'); xlabel('w'); ylabel('h(n)'); grid("on");


[20/11/25, 3:42:55 AM] Cocomelon: 



clc; clear; close all
N = input('Enter the number of elements: ')
wc1 = input('Enter the lower cutoff frequency: ')
wc2 = input('Enter the higher cutoff frequency: ')

a = (N - 1) / 2;
n = -a:a;
mid = a + 1;

h = (sin(wc2*n) - sin(wc1*n)) ./ (pi*n);
h(mid) = (wc2 - wc1) / pi;

disp('Impulse Response h(n):');
disp(h);

h_r  = h .* 1
h_hn = h .* (0.5 + 0.5 * cos((pi*n) ./ a))
h_hm = h .* (0.54 + 0.46 * cos((pi*n) ./ a))
h_bl = h .* (0.42 + 0.5 * cos((pi*n) ./ a) + 0.08 * cos((2*pi*n) ./ a))

[Hr, w] = freqz(h_r, 1, 1024);
[Hn, ~] = freqz(h_hn, 1, 1024);
[Hm, ~] = freqz(h_hm, 1, 1024);
[Hb, ~] = freqz(h_bl, 1, 1024);

figure('Name','Band Pass FIR Filter');
subplot(2,2,1); plot(w, abs(Hr)); title('Rectangular Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,2); plot(w, abs(Hn)); title('Hanning Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,3); plot(w, abs(Hm)); title('Hamming Window'); xlabel('w'); ylabel('h(n)'); grid("on");
subplot(2,2,4); plot(w, abs(Hb)); title('Blackman Window'); xlabel('w'); ylabel('h(n)'); grid("on");

