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

N = input('N: ');
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

N = input('N: ');
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

N = input('N: ');
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

