clc; clear; close all;

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