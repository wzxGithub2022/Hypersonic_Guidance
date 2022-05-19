%% 使用三阶切比雪夫多项式在区间[-pi/2,pi/2]上逼近sin
% f(x) = c0/2 + sum(ck*Tk(x)); k = 1,2,...
% ck = 2/pi * integral(f(cos(theta)) * cos(k*theta)) d_theta; [0,pi]
% x = pi/2 * t; t belongs to [-1,+1]
% t = x/(pi/2);

clear;clc;
syms theta
c0_half = (1/pi) * int( sin( pi/2 * cos(theta) ),0,pi);
c0h_double = double(c0_half);

for k = 1:3

    c(k) = (2/pi) * int(cos(k*theta) * sin( pi/2 * cos(theta) ) ,0,pi);
    c_double(k) = double(c(k));
    
end
% 1.133648177811748 -0.138071776587192
% T0_x = 1;
% T1_x = @(x) x;
% T2_x = @(x) 2*x^2 - 1;
% T3_x = @(x) 4*x^3 - 3*x;
% T4_x = @(x) 8*x^4 - 8*x^2 + 1;

ft = @(t) 1.133648177811748 * t - 0.138071776587192 * (4*t.^3 - 3*t);

t = linspace(-1,1,100);
sin_cheby = ft(t);

x = pi/2 * t;
y = sin(x);

figure(1)
subplot(121),plot(x,sin_cheby,x,y);
legend('ft','original');
title('切比雪夫三阶多项式逼近效果');

%%
syms t x

t = x/(pi/2);
ft = 1.133648177811748 * t - 0.138071776587192 * (4*t.^3 - 3*t);
fx = expand(ft);
disp(fx);
% (27883830063710443*x)/(9007199254740992*pi) - (2487280006353841*x^3)/(562949953421312*pi^3)

a0 = 0;
a1 = (27883830063710443)/(9007199254740992*pi);     % 0.985400513847416
a2 = 0;
a3 = - (2487280006353841)/(562949953421312*pi^3);   % -0.142496853019355

x = linspace(-pi/2,pi/2,100);
y = sin(x);

gx = a0 + a1*x + a2*x.^2 + a3*x.^3;

figure(1)
subplot(122),plot(x,gx,x,y);
legend('gx','original');
title('切比雪夫三阶多项式逼近效果');

%%
x = linspace(-pi,pi,100);
y = sin(x);
gx = a0 + a1*x + a2*x.^2 + a3*x.^3;

tl = -pi/2 * ones(100,1);
tu = pi/2 * ones(100,1);
yl = linspace(-2,2,100);
yu = linspace(-2,2,100);

figure(2)
plot(x,gx,'-*r',x,y,'-ob');
hold on
plot(tl,yl,'--k',tu,yu,'--k','LineWidth',2);
axis([-pi pi -2 2]);
grid on
legend('多项式','正弦');
title('切比雪夫三阶多项式逼近效果');



