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

% ===== NEW: Continuous Trigonometric Signals =====
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