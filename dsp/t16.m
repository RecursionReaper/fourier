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