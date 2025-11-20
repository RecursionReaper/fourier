clc; clear; close all;

x = input('Enter sequence x: ')

n = input('Enter index vector n: ')

a = input('Enter scaling factor a: ')

m = input('Enter shift value m: ')



n1 = n - m;
figure(1);
stem(n1 / a, x .* (mod(n1, a) == 0), 'filled'); 
title('x[a n + m] (Shift + Scale)'); 
xlabel('n'); ylabel('Amplitude'); grid on;