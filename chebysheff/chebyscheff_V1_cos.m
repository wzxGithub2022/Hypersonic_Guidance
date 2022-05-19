%% 使用切比雪夫多项式在区间[-pi/2,pi/2]上逼近cos
% f(x) = c0/2 + sum(ck*Tk(x)); k = 1,2,...
% ck = 2/pi * integral(f(cos(theta)) * cos(k*theta)) d_theta; [0,pi]
% x = pi/2 * t; t belongs to [-1,+1]
% t = x/(pi/2);

clear;clc;
syms theta
c0_half = (1/pi) * int( cos( pi/2 * cos(theta) ),0,pi);
c0h_double = double(c0_half);

for k = 1:4

    c(k) = (2/pi) * int(cos(k*theta) * cos( pi/2 * cos(theta) ) ,0,pi);
    c_double(k) = double(c(k));
    
end
% 0.472001215768235 -0.499403258270407 0.027992079617548
% T0_x = 1;
% T1_x = @(x) x;
% T2_x = @(x) 2*x^2 - 1;
% T3_x = @(x) 4*x^3 - 3*x;
% T4_x = @(x) 8*x^4 - 8*x^2 + 1;

ft = @(t) 0.472001215768235 - 0.499403258270407 * (2*t.^2 - 1) + 0.027992079617548 * (8*t.^4 - 8*t.^2 + 1);

t = linspace(-1,1,100);
cos_cheby = ft(t);

x = pi/2 * t;
y = cos(x);

figure(1)
subplot(121),plot(x,cos_cheby,x,y);
legend('ft','original');
title('切比雪夫多项式逼近效果');

%%
syms t x

t = x/(pi/2);
ft = 0.472001215768235 - 0.499403258270407 * (2*t.^2 - 1) + 0.027992079617548 * (8*t.^4 - 8*t.^2 + 1);
fx = expand(ft);
disp(fx);
% (8068167637434523*x^4)/(2251799813685248*pi^4) - (44053964883101987*x^2)/(9007199254740992*pi^2) + 288056444585048011/288230376151711744

a0 = 288056444585048011/288230376151711744;         % 0.999396553656190
a1 = 0;     
a2 = - (44053964883101987)/(9007199254740992*pi^2); % -0.495559134405119
a3 = 0;   
a4 = (8068167637434523)/(2251799813685248*pi^4);    % 0.036782872656059

x = linspace(-pi/2,pi/2,100);
y = cos(x);

gx = a0 + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4;

figure(1)
subplot(122),plot(x,gx,x,y);
legend('gx','original');
title('切比雪夫多项式逼近效果');

%%
x = linspace(-pi,pi,100);
y = cos(x);
gx = a0 + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4;

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
legend('多项式','余弦');
title('切比雪夫三阶多项式逼近效果');
