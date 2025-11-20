clc; clear; close all;

Ac = input('Ac: ');
t0 = input('t0: ');
Ad = input('Ad: ');
n0 = input('n0: ');

t = -10:0.1:10;
imp_c = (t==0);
stp_c = (t>=0);
rmp_c = t.*(t>=0);
exp_c = exp(t);

figure;
subplot(2,2,1); plot(t,imp_c); title('impulse(t)');
subplot(2,2,2); plot(t,stp_c); title('u(t)');
subplot(2,2,3); plot(t,rmp_c); title('r(t)');
subplot(2,2,4); plot(t,exp_c); title('e^t');

tt = t - t0;

figure;
subplot(2,2,1); plot(t,Ac*(tt==0)); title('impulse(t-t0)');
subplot(2,2,2); plot(t,Ac*(tt>=0)); title('u(t-t0)');
subplot(2,2,3); plot(t,Ac*(tt.*(tt>=0))); title('r(t-t0)');
subplot(2,2,4); plot(t,Ac*exp(tt)); title('e^{t-t0}');

w = pi/3;
y1 = sin(w*t);
y2 = cos(w*t);
y3 = sin(w*t)./(w*t); y3(t==0)=1;

figure;
subplot(3,1,1); plot(t,y1); title('sin(wt)');
subplot(3,1,2); plot(t,y2); title('cos(wt)');
subplot(3,1,3); plot(t,y3); title('sinc(wt)');

yc1 = Ac*sin(w*tt);
yc2 = Ac*cos(w*tt);
yc3 = Ac*(sin(w*tt)./(w*tt)); yc3(tt==0)=Ac;

figure;
subplot(3,1,1); plot(t,yc1); title('sin(w(t-t0))');
subplot(3,1,2); plot(t,yc2); title('cos(w(t-t0))');
subplot(3,1,3); plot(t,yc3); title('sinc(w(t-t0))');

n = -10:10;
imp = (n==0);
stp = (n>=0);
rmp = n.*(n>=0);
expd = exp(n);

figure;
subplot(2,2,1); stem(n,imp,'filled'); title('impulse[n]');
subplot(2,2,2); stem(n,stp,'filled'); title('u[n]');
subplot(2,2,3); stem(n,rmp,'filled'); title('r[n]');
subplot(2,2,4); stem(n,expd,'filled'); title('e^n');

nn = n - n0;

figure;
subplot(2,2,1); stem(n,Ad*(nn==0),'filled'); title('impulse[n-n0]');
subplot(2,2,2); stem(n,Ad*(nn>=0),'filled'); title('u[n-n0]');
subplot(2,2,3); stem(n,Ad*(nn.*(nn>=0)),'filled'); title('r[n-n0]');
subplot(2,2,4); stem(n,Ad*exp(nn),'filled'); title('e^{n-n0}');

y1d = sin(w*n);
y2d = cos(w*n);
y3d = sin(w*n)./(w*n); y3d(n==0)=1;

figure;
subplot(3,1,1); stem(n,y1d,'filled'); title('sin(wn)');
subplot(3,1,2); stem(n,y2d,'filled'); title('cos(wn)');
subplot(3,1,3); stem(n,y3d,'filled'); title('sinc(wn)');

y1dc = Ad*sin(w*nn);
y2dc = Ad*cos(w*nn);
y3dc = Ad*(sin(w*nn)./(w*nn)); y3dc(nn==0)=Ad;

figure;
subplot(3,1,1); stem(n,y1dc,'filled'); title('sin(w(n-n0))');
subplot(3,1,2); stem(n,y2dc,'filled'); title('cos(w(n-n0))');
subplot(3,1,3); stem(n,y3dc,'filled'); title('sinc(w(n-n0))');